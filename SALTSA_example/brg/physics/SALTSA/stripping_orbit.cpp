/**********************************************************************\
  @file stripping_orbit.cpp

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

#include <cstdlib>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <utility>

#include "brg/global.h"

#include "brg/file_functions.h"

#include "brg/math/calculus/integrate.hpp"
#include "brg/math/interpolator/interpolator.h"
#include "brg/math/interpolator/interpolator_derivative.h"
#include "brg/math/solvers/solvers.hpp"

#include "brg/physics/astro.h"
#include "brg/physics/SALTSA/gabdt.h"
#include "brg/physics/SALTSA/solve_rt_functors.h"
#include "brg/physics/phase.hpp"
#include "brg/physics/units/unit_obj.h"

#include "stripping_orbit.h"

using namespace std;

/** Static member initialisations **/
#if (1)

// stripping_orbit static member initialisations
#if (1)
// Default integration parameters
#if(1)
// Default number of steps for which stripping is calculated
size_t brgastro::stripping_orbit::_default_spline_resolution_ = 100;

// Default interpolation type
brgastro::stripping_orbit::allowed_interpolation_type brgastro::stripping_orbit::_default_interpolation_type_ = UNSET;

// Variable step length tweaking: Time step length is proportional to (v_0/v)^(step_length_power)
// This gives smaller steps when the satellite is moving faster.
// If you want to turn off adaptive step size, set step_length_power to 0
// Alternatively, set step_length_power to 1 for even steps in position
BRG_VELOCITY brgastro::stripping_orbit::_default_v_0_ = 400 * unitconv::kmpstomps; // 400 km/s
BRG_DISTANCE brgastro::stripping_orbit::_default_r_0_ = 400 * unitconv::kpctom; // 400 kpc
double brgastro::stripping_orbit::_default_step_length_power_ = 3.0;
double brgastro::stripping_orbit::_default_step_factor_max_ = 10; // Maximum allowed value of (v_0/v)^(step_length_power)
double brgastro::stripping_orbit::_default_step_factor_min_ = 0.001; // Minimum allowed value of (v_0/v)^(step_length_power)
#endif

// Default tuning values

#if(1)
// Tuning parameters, for how strong stripping and shocking are and when shocking is active
double brgastro::stripping_orbit::_default_tidal_stripping_amplification_ = 0.666; // Tuned
double brgastro::stripping_orbit::_default_tidal_stripping_deceleration_ = 0.; // Tuned
double brgastro::stripping_orbit::_default_tidal_stripping_radialness_ = 0.; // Tuned
double brgastro::stripping_orbit::_default_tidal_shocking_amplification_ = 3.0; // Tuned
double brgastro::stripping_orbit::_default_tidal_shocking_persistance_ = 1.0; // How long shocking is active for
double brgastro::stripping_orbit::_default_tidal_shocking_power_ = -1.5; // Affects interplay of stripping and satellite halo profile
#endif

#endif

#endif

/** Class method implementations **/
#if (1)

// brgastro::stripping_orbit class method implementations
#if (1)

void brgastro::stripping_orbit::_pass_parameters_to_segment(
		brgastro::stripping_orbit_segment & segment,
		brgastro::density_profile *segment_init_satellite,
		brgastro::density_profile *segment_init_host,
		size_t resolution) const
{
	if(segment_init_satellite==NULL)
	{
		segment.set_init_satellite(
				_init_satellite_ptr_ );
	}
	else
	{
		segment.set_init_satellite(
				segment_init_satellite );
	}

	if(segment_init_host==NULL)
	{
		segment.set_init_host( _init_host_ptr_ );
	}
	else
	{
		segment.set_init_host( segment_init_host );
	}

	if(resolution==0)
	{
		segment.set_resolution( _base_resolution_);
	}
	else
	{
		segment.set_resolution( resolution );
	}

	segment.set_interpolation_type(_interpolation_type_);
	segment.set_v_0(_v_0_ );
	segment.set_r_0(_r_0_ );
	segment.set_step_length_power(_step_length_power_ );
	segment.set_step_factor_min(_step_factor_min_ );
	segment.set_step_factor_max(_step_factor_max_ );
	segment.set_tidal_stripping_amplification(_tidal_stripping_amplification_ );
	segment.set_tidal_stripping_deceleration(_tidal_stripping_deceleration_ );
	segment.set_tidal_stripping_radialness(_tidal_stripping_radialness_ );
	segment.set_tidal_shocking_amplification(_tidal_shocking_amplification_ );
	segment.set_tidal_shocking_persistance(_tidal_shocking_persistance_ );
	segment.set_tidal_shocking_power(_tidal_shocking_power_ );
	segment.set_record_full_data(_record_full_data_ );

	return;
}

const std::vector< brgastro::stripping_orbit_segment >::iterator brgastro::stripping_orbit::_final_good_segment() const
{
	if(!_calculated_)
	{
		_final_good_segment_ = _orbit_segments_.end();
		calc(); // It might not work, but running it will at least give us the most sensible result
	}
	return _final_good_segment_;
}

// Swap functions
void brgastro::stripping_orbit::swap(stripping_orbit &other)
{
	using std::swap;

	// Integration parameters
#if(1)
	swap(_base_resolution_, other._base_resolution_);
	swap(_interpolation_type_, other._interpolation_type_);
	swap(_v_0_, other._v_0_);
	swap(_r_0_, other._r_0_);
	swap(_step_length_power_, other._step_length_power_);
	swap(_step_factor_max_, other._step_factor_max_);
	swap(_step_factor_min_, other._step_factor_min_);
#endif

	// Tuning values
#if(1)
	swap(_tidal_stripping_amplification_, other._tidal_stripping_amplification_);
	swap(_tidal_stripping_deceleration_, other._tidal_stripping_deceleration_);
	swap(_tidal_stripping_radialness_, other._tidal_stripping_radialness_);
	swap(_tidal_shocking_amplification_, other._tidal_shocking_amplification_);
	swap(_tidal_shocking_persistance_, other._tidal_shocking_persistance_);
	swap(_tidal_shocking_power_, other._tidal_shocking_power_);
#endif

	swap(_satellite_parameter_unitconvs_,
			other._satellite_parameter_unitconvs_);
	swap(_satellite_output_parameters_,
			other._satellite_output_parameters_);
	swap(_host_parameter_unitconvs_,
			other._host_parameter_unitconvs_);
	swap(_host_output_parameters_, other._host_output_parameters_);
	swap(_record_full_data_, other._record_full_data_);
	swap(_num_segments_, other._num_segments_);

	swap(_t_min_natural_value_, other._t_min_natural_value_);
	swap(_t_max_natural_value_, other._t_max_natural_value_);
	swap(_t_min_override_value_, other._t_min_override_value_);
	swap(_t_max_override_value_, other._t_max_override_value_);
	swap(_override_t_min_, other._override_t_min_);
	swap(_override_t_max_, other._override_t_max_);

	swap(_final_frac_m_ret_list_, other._final_frac_m_ret_list_);
	swap(_final_frac_m_vir_ret_list_, other._final_frac_m_vir_ret_list_);

	swap(_x_points_, other._x_points_);
	swap(_y_points_, other._y_points_);
	swap(_z_points_, other._z_points_);
	swap(_test_mass_points_,
			other._test_mass_points_);
	_test_mass_error_interpolator_.swap(
			other._test_mass_error_interpolator_);
	_test_mass_interpolator_.swap(
			other._test_mass_interpolator_);
	swap(_vx_points_, other._vx_points_);
	swap(_vy_points_, other._vy_points_);
	swap(_vz_points_, other._vz_points_);
	swap(_vx_unknown_points_,
			other._vx_unknown_points_);
	swap(_vy_unknown_points_,
			other._vy_unknown_points_);
	swap(_vz_unknown_points_,
			other._vz_unknown_points_);
	swap(_host_parameter_points_,
			other._host_parameter_points_);
	_m_ret_interpolator_.swap(
			other._m_ret_interpolator_);
	_m_vir_ret_interpolator_.swap(
			other._m_vir_ret_interpolator_);

	swap(_t_points_, other._t_points_);
	swap(_host_param_t_points_, other._host_param_t_points_);

	swap(_discontinuity_times_, other._discontinuity_times_);
	swap(_cleaned_discontinuity_times_,
			other._cleaned_discontinuity_times_);

	swap(_init_host_ptr_, other._init_host_ptr_);
	swap(_init_satellite_ptr_, other._init_satellite_ptr_);

	swap(_host_is_evolving_, other._host_is_evolving_);
	swap(_calculated_, other._calculated_);
	swap(_bad_result_, other._bad_result_);
	swap(_host_loaded_, other._host_loaded_);
	swap(_satellite_loaded_, other._satellite_loaded_);
	swap(_using_private_init_host_,
			other._using_private_init_host_);
	swap(_using_private_init_satellite_,
			other._using_private_init_satellite_);
	swap(_private_tNFW_init_host_, other._private_tNFW_init_host_);
	swap(_private_tNFW_init_satellite_,
			other._private_tNFW_init_satellite_);

	swap(_orbit_segments_, other._orbit_segments_);
	swap(_final_good_segment_,other._final_good_segment_);

	swap(_likely_disrupted_, other._likely_disrupted_);

	// It's possible the addresses of _private_tNFW_init_host_ and _private_tNFW_init_satellite_
	// changed, so correct the pointers to them if they're in use

	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( other._using_private_init_host_ )
	{
		other._init_host_ptr_ = &(other._private_tNFW_init_host_);
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
	if ( other._using_private_init_satellite_ )
	{
		other._init_satellite_ptr_ = &(other._private_tNFW_init_satellite_);
	}
}
namespace std
{
	template <>
	void swap(brgastro::stripping_orbit &same,
			brgastro::stripping_orbit &other)
	{
		same.swap(other);
	}
}

brgastro::stripping_orbit::stripping_orbit()
{
	clear();
}

brgastro::stripping_orbit::stripping_orbit(
		const stripping_orbit &other )
{
	// Integration parameters
#if(1)
	_base_resolution_ = other._base_resolution_;
	_interpolation_type_ = other._interpolation_type_;
	_v_0_ = other._v_0_;
	_r_0_ = other._r_0_;
	_step_length_power_ = other._step_length_power_;
	_step_factor_max_ = other._step_factor_max_;
	_step_factor_min_ = other._step_factor_min_;
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = other._tidal_stripping_amplification_;
	_tidal_stripping_deceleration_ = other._tidal_stripping_deceleration_;
	_tidal_stripping_radialness_ = other._tidal_stripping_radialness_;
	_tidal_shocking_amplification_ = other._tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = other._tidal_shocking_persistance_;
	_tidal_shocking_power_ = other._tidal_shocking_power_;
#endif

	_satellite_parameter_unitconvs_ =
			other._satellite_parameter_unitconvs_;
	_satellite_output_parameters_ =
			other._satellite_output_parameters_;
	_host_parameter_unitconvs_ =
			other._host_parameter_unitconvs_;
	_host_output_parameters_ = other._host_output_parameters_;
	_record_full_data_ = other._record_full_data_;
	_num_segments_ = other._num_segments_;

	_t_min_natural_value_ = other._t_min_natural_value_;
	_t_max_natural_value_ = other._t_max_natural_value_;
	_t_min_override_value_ = other._t_min_override_value_;
	_t_max_override_value_ = other._t_max_override_value_;
	_override_t_min_ = other._override_t_min_;
	_override_t_max_ = other._override_t_max_;

	_final_frac_m_ret_list_ = other._final_frac_m_ret_list_;
	_final_frac_m_vir_ret_list_ = other._final_frac_m_vir_ret_list_;

	_x_points_ = other._x_points_;
	_y_points_ = other._y_points_;
	_z_points_ = other._z_points_;
	_test_mass_points_ =
			other._test_mass_points_;
	_test_mass_error_interpolator_ =
			other._test_mass_error_interpolator_;
	_test_mass_interpolator_ =
			other._test_mass_interpolator_;
	_vx_points_ = other._vx_points_;
	_vy_points_ = other._vy_points_;
	_vz_points_ = other._vz_points_;
	_vx_unknown_points_ =
			other._vx_unknown_points_;
	_vy_unknown_points_ =
			other._vy_unknown_points_;
	_vz_unknown_points_ =
			other._vz_unknown_points_;
	_host_parameter_points_ =
			other._host_parameter_points_;
	_m_ret_interpolator_ =
			other._m_ret_interpolator_;
	_m_vir_ret_interpolator_ =
			other._m_vir_ret_interpolator_;

	_t_points_ = other._t_points_;
	_host_param_t_points_ = other._host_param_t_points_;

	_discontinuity_times_ = other._discontinuity_times_;
	_cleaned_discontinuity_times_ =
			other._cleaned_discontinuity_times_;

	_init_host_ptr_ = other._init_host_ptr_;
	_init_satellite_ptr_ = other._init_satellite_ptr_;
	_host_is_evolving_ = other._host_is_evolving_;
	_calculated_ = other._calculated_;
	_bad_result_ = other._bad_result_;
	_host_loaded_ = other._host_loaded_;
	_satellite_loaded_ = other._satellite_loaded_;
	_using_private_init_host_ =
			other._using_private_init_host_;
	_using_private_init_satellite_ =
			other._using_private_init_satellite_;
	_private_tNFW_init_host_ = other._private_tNFW_init_host_;
	_private_tNFW_init_satellite_ =
			other._private_tNFW_init_satellite_;

	_orbit_segments_ = other._orbit_segments_;

	_likely_disrupted_ = other._likely_disrupted_;

	if( (other._final_good_segment_ == other._orbit_segments_.end()) || ( _orbit_segments_.empty() ) )
	{
		_final_good_segment_ = _orbit_segments_.end();
	}
	else
	{
		// We'll have to find the segment
		std::vector< brgastro::stripping_orbit_segment >::iterator test_iterator_old = other._orbit_segments_.begin();
		std::vector< brgastro::stripping_orbit_segment >::iterator test_iterator_new = _orbit_segments_.begin();
		bool found = false;
		while(test_iterator_old != other._orbit_segments_.end())
		{
			if(test_iterator_old == other._final_good_segment_)
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

brgastro::stripping_orbit & brgastro::stripping_orbit::operator=(
		stripping_orbit other )
{
	swap(other);
	return *this;
}

brgastro::stripping_orbit::stripping_orbit( const density_profile *init_init_host,
		const density_profile *init_init_satellite, const int init_resolution )
{
	clear();
	_init_host_ptr_ = init_init_host;
	_init_satellite_ptr_ = init_init_satellite;
	_base_resolution_ = init_resolution;
	_host_loaded_ = true;
	_satellite_loaded_ = true;

}

brgastro::stripping_orbit::~stripping_orbit()
{
}

brgastro::stripping_orbit *brgastro::stripping_orbit::stripping_orbit_clone()
{
	return new stripping_orbit( *this );
}

void brgastro::stripping_orbit::clear()
{
	clear_calcs();
	clear_points();

	// Integration parameters
#if(1)
	_base_resolution_ = _default_spline_resolution_;
	_interpolation_type_ = _default_interpolation_type_;
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
	_tidal_stripping_radialness_ = _default_tidal_stripping_radialness_;
	_tidal_shocking_amplification_ = _default_tidal_shocking_amplification_;
	_tidal_shocking_persistance_ = _default_tidal_shocking_persistance_;
	_tidal_shocking_power_ = _default_tidal_shocking_power_;
#endif

	_satellite_parameter_unitconvs_.clear();
	_satellite_output_parameters_.clear();
	_host_parameter_unitconvs_.clear();
	_host_output_parameters_.clear();

	_t_min_override_value_ = DBL_MAX;
	_t_max_override_value_ = ( -DBL_MAX );
	_override_t_min_ = false;
	_override_t_max_ = false;
	_record_full_data_ = false;

	_init_host_ptr_ = 0;
	_init_satellite_ptr_ = 0;

	_host_loaded_ = false;
	_satellite_loaded_ = false;
	_using_private_init_host_ = false;
	_using_private_init_satellite_ = false;
	_host_is_evolving_ = false;
}

void brgastro::stripping_orbit::clear_calcs() const
{
	_orbit_segments_.clear();
	_num_segments_ = 0;
	_final_good_segment_ = _orbit_segments_.end();
	_final_frac_m_ret_list_.clear();
	_final_frac_m_vir_ret_list_.clear();
	_m_ret_interpolator_.clear();
	_m_vir_ret_interpolator_.clear();

	_calculated_ = false;
	_bad_result_ = false;
	_likely_disrupted_ = false;
}

void brgastro::stripping_orbit::add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
		const double test_mass, const double test_mass_error )
{
	// Check if there's already a point with this t value
	for(size_t i=0; i<_t_points_.size(); i++)
	{
		if(t==_t_points_[i])
			throw std::runtime_error("Attempt to add duplicate t value to stripping_orbit.");
	}
	force_add_point(x, y, z, t, test_mass, test_mass_error );
}

void brgastro::stripping_orbit::force_add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
		const double test_mass, const double test_mass_error )
{
	_calculated_ = false;
	try
	{
		_x_points_.push_back( std::pair< double, double >( t, x ) );
		_y_points_.push_back( std::pair< double, double >( t, y ) );
		_z_points_.push_back( std::pair< double, double >( t, z ) );
		_vx_unknown_points_.push_back( t );
		_vy_unknown_points_.push_back( t );
		_vz_unknown_points_.push_back( t );
		_test_mass_points_.push_back(
				std::pair< double, double >( t, test_mass ) );
		_test_mass_interpolator_.add_point(t, test_mass);
		_test_mass_error_interpolator_.add_point(t, test_mass_error);
		_t_points_.push_back(t);
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n";
		throw;
	}
}

void brgastro::stripping_orbit::add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
		CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
		const double test_mass, const double test_mass_error )
{
	// Check if there's already a point with this t value
	for(size_t i=0; i<_t_points_.size(); i++)
	{
		if(t==_t_points_[i])
			throw std::runtime_error("Attempt to add duplicate t value to stripping_orbit.");
	}
	return force_add_point(x, y, z, vx, vy, vz, t, test_mass, test_mass_error );
}

void brgastro::stripping_orbit::force_add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
		CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
		const double test_mass, const double test_mass_error )
{
	_calculated_ = false;
	try
	{
		_x_points_.push_back( std::pair< double, double >( t, x ) );
		_y_points_.push_back( std::pair< double, double >( t, y ) );
		_z_points_.push_back( std::pair< double, double >( t, z ) );
		_vx_points_.push_back( std::pair< double, double >( t, vx ) );
		_vy_points_.push_back( std::pair< double, double >( t, vy ) );
		_vz_points_.push_back( std::pair< double, double >( t, vz ) );
		_test_mass_points_.push_back(
				std::pair< double, double >( t, test_mass ) );
		_test_mass_interpolator_.add_point( t, test_mass );
		_test_mass_error_interpolator_.add_point( t, test_mass_error );
		_t_points_.push_back(t);
		if ( t < _t_min_natural_value_ )
			_t_min_natural_value_ = t;
		if ( t > _t_max_natural_value_ )
			_t_max_natural_value_ = t;
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_point().\n";
		throw;
	}
}

void brgastro::stripping_orbit::add_host_parameter_point(
		const std::vector< BRG_UNITS > &parameters, CONST_BRG_TIME_REF t,
		const bool silent )
{
	// Check if there's already a point with this t value
	for(size_t i=0; i<_host_param_t_points_.size(); i++)
	{
		if(t==_host_param_t_points_[i])
			throw std::runtime_error("Attempt to add duplicate t value to stripping_orbit.");
	}
	force_add_host_parameter_point( parameters, t, silent );
}

void brgastro::stripping_orbit::force_add_host_parameter_point(
		const std::vector< BRG_UNITS > &parameters, CONST_BRG_TIME_REF t,
		const bool silent )
{
	// Check num_parameters matches vector size
	if ( _init_host_ptr_->num_parameters() != parameters.size() )
	{
		throw std::logic_error("ERROR: num_parameters must == parameters.size() in stripping_orbit::add_host_parameter_point.\n");
	}

	_host_parameter_points_.push_back(
			std::pair< double, std::vector< BRG_UNITS > >( t, parameters ) );
	_host_param_t_points_.push_back(t);

	if ( _host_parameter_points_.size() >= 2 )
		_host_is_evolving_ = true;

	_calculated_ = false;
}

void brgastro::stripping_orbit::add_discontinuity_time(
		CONST_BRG_TIME_REF t )
{
	try
	{
		_discontinuity_times_.push_back( t );
	}
	catch (const std::exception &e)
	{
		std::cerr << "ERROR: Exception in stripping_orbit::add_discontinuity_time().\n";
		throw;
	}
}

void brgastro::stripping_orbit::clear_points()
{
	_x_points_.clear();
	_y_points_.clear();
	_z_points_.clear();

	_vx_points_.clear();
	_vy_points_.clear();
	_vz_points_.clear();

	_vx_unknown_points_.clear();
	_vy_unknown_points_.clear();
	_vz_unknown_points_.clear();

	_test_mass_points_.clear();
	_test_mass_interpolator_.clear();
	_test_mass_error_interpolator_.clear();

	_t_points_.clear();
	_host_param_t_points_.clear();

	_t_min_natural_value_ = DBL_MAX;
	_t_max_natural_value_ = ( -DBL_MAX );

	_calculated_ = false;
}
void brgastro::stripping_orbit::clear_discontinuity_times()
{
	_discontinuity_times_.clear();
	_calculated_ = false;
}
void brgastro::stripping_orbit::clear_host_parameter_points()
{
	_host_parameter_points_.clear();
	_host_param_t_points_.clear();
	_calculated_ = false;
}

void brgastro::stripping_orbit::set_init_host(
		const density_profile *new_init_host )
{
	_init_host_ptr_ = new_init_host;
	_host_loaded_ = true;
	_using_private_init_host_ = false;
	_calculated_ = false;
}

void brgastro::stripping_orbit::set_init_satellite(
		const density_profile *new_init_satellite )
{
	_init_satellite_ptr_ = new_init_satellite;
	_satellite_loaded_ = true;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
}

void brgastro::stripping_orbit::set_tNFW_init_satellite(
		CONST_BRG_MASS_REF new_init_mvir0, const double z,
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
}

void brgastro::stripping_orbit::set_tNFW_init_host(
		CONST_BRG_MASS_REF new_init_mvir0, const double z,
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
}

void brgastro::stripping_orbit::clear_init_satellite()
{
	_init_satellite_ptr_ = 0;
	_satellite_loaded_ = false;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
}
void brgastro::stripping_orbit::clear_init_host()
{
	_init_host_ptr_ = 0;
	_host_loaded_ = false;
	_using_private_init_host_ = false;
	_calculated_ = false;
}

// Setting default integration parameters
#if(1)
void brgastro::stripping_orbit::set_default_resolution( const size_t new_default_resolution)
{
	if ( new_default_resolution < 2 )
	{
		throw std::logic_error("ERROR: Attempt to set default resolution to value below minimum of 2.\n");
	}
	_default_spline_resolution_ = new_default_resolution;
}
void brgastro::stripping_orbit::set_default_resolution( const size_t new_default_resolution,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_resolution == _default_spline_resolution_ )
		return;

	if ( new_default_resolution < 2 )
	{
		throw std::logic_error("WARNING: Attempt to set default resolution to value below minimum of 2.\n");
	}
	_default_spline_resolution_ = new_default_resolution;

	if(override_current) reset_resolution();
}
void brgastro::stripping_orbit::set_default_interpolation_type(
		const allowed_interpolation_type new_default_interpolation_type)
{
	_default_interpolation_type_ = new_default_interpolation_type;
}
void brgastro::stripping_orbit::set_default_interpolation_type(
		const allowed_interpolation_type new_default_interpolation_type,
		const bool override_current,
		const bool silent)
{
	// Check if anything is actually changing here
	if ( new_default_interpolation_type == _default_interpolation_type_ )
		return;
	_default_interpolation_type_ = new_default_interpolation_type;

	if(override_current) reset_interpolation_type();
}
void brgastro::stripping_orbit::set_default_v_0( const double new_default_v_0)
{

	if ( new_default_v_0 <= 0 )
	{
		throw std::logic_error("WARNING: Attempt to set default v_0 to value at or below minimum of 0.\n");
	}
	_default_v_0_ = new_default_v_0;
}
void brgastro::stripping_orbit::set_default_v_0( const double new_default_v_0,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_v_0 == _default_v_0_ )
		return;

	if ( new_default_v_0 <= 0 )
	{
		throw std::logic_error("WARNING: Attempt to set default v_0 to value at or below minimum of 0.\n");
	}
	_default_v_0_ = new_default_v_0;

	if(override_current) reset_v_0();
}
void brgastro::stripping_orbit::set_default_r_0( const double new_default_r_0)
{

	if ( new_default_r_0 <= 0 )
	{
		throw std::logic_error("Attempt to set default r_0 to value at or below minimum of 0.");
	}
	_default_r_0_ = new_default_r_0;
}
void brgastro::stripping_orbit::set_default_r_0( const double new_default_r_0,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_r_0 == _default_r_0_ )
		return;

	if ( new_default_r_0 <= 0 )
	{
		throw std::logic_error("Attempt to set default r_0 to value at or below minimum of 0.");
	}
	_default_r_0_ = new_default_r_0;

	if(override_current) reset_r_0();
}
void brgastro::stripping_orbit::set_default_step_length_power(
		const double new_default_step_length_power)
{
	_default_step_length_power_ = new_default_step_length_power;
}
void brgastro::stripping_orbit::set_default_step_length_power(
		const double new_default_step_length_power,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_length_power == _default_step_length_power_ )
		return;
	_default_step_length_power_ = new_default_step_length_power;

	if(override_current) reset_step_length_power();
}
void brgastro::stripping_orbit::set_default_step_factor_max(
		const double new_default_step_factor_max)
{

	if ( new_default_step_factor_max < 1 )
	{
		throw std::logic_error("Attempt to set default step_factor_max to value below minimum of 1.");
	}
	_default_step_factor_max_ = new_default_step_factor_max;
}
void brgastro::stripping_orbit::set_default_step_factor_max(
		const double new_default_step_factor_max,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_factor_max == _default_step_factor_max_ )
		return;

	if ( new_default_step_factor_max < 1 )
	{
		throw std::logic_error("Attempt to set default step_factor_max to value below minimum of 1.");
	}
	_default_step_factor_max_ = new_default_step_factor_max;

	if(override_current) reset_step_factor_max();
}
void brgastro::stripping_orbit::set_default_step_factor_min(
		const double new_default_step_factor_min)
{
	_default_step_factor_min_ = new_default_step_factor_min;
}
void brgastro::stripping_orbit::set_default_step_factor_min(
		const double new_default_step_factor_min,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_step_factor_min == _default_step_factor_min_ )
		return;

	if ( new_default_step_factor_min > 1 )
	{
		throw std::logic_error("Attempt to set default step_factor_min to value above minimum of 1.");
	}
	_default_step_factor_min_ = new_default_step_factor_min;

	if(override_current) reset_step_factor_min();
}
#endif

// Setting default tuning parameters
#if(1)

void brgastro::stripping_orbit::set_default_tidal_stripping_amplification(
		const double new_default_tidal_stripping_amplification)
{
	_default_tidal_stripping_amplification_ = new_default_tidal_stripping_amplification;
}

void brgastro::stripping_orbit::set_default_tidal_stripping_amplification(
		const double new_default_tidal_stripping_amplification,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_stripping_amplification == _default_tidal_stripping_amplification_ )
		return;

	if ( new_default_tidal_stripping_amplification < 0 )
	{
		throw std::logic_error("Attempt to set default tidal_stripping_amplification to value below minimum of 0.");
	}
	_default_tidal_stripping_amplification_ = new_default_tidal_stripping_amplification;

	if(override_current)
		reset_tidal_stripping_amplification();
}
void brgastro::stripping_orbit::set_default_tidal_stripping_deceleration(
		const double new_default_tidal_stripping_deceleration)
{
	_default_tidal_stripping_deceleration_ = new_default_tidal_stripping_deceleration;
}
void brgastro::stripping_orbit::set_default_tidal_stripping_deceleration(
		const double new_default_tidal_stripping_deceleration,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_stripping_deceleration == _default_tidal_stripping_deceleration_ )
		return;
	_default_tidal_stripping_deceleration_ = new_default_tidal_stripping_deceleration;

	if(override_current)
			reset_tidal_stripping_deceleration();
}
void brgastro::stripping_orbit::set_default_tidal_stripping_radialness(
		const double new_default_tidal_stripping_radialness)
{
	_default_tidal_stripping_radialness_ = new_default_tidal_stripping_radialness;
}
void brgastro::stripping_orbit::set_default_tidal_stripping_radialness(
		const double new_default_tidal_stripping_radialness,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_stripping_radialness == _default_tidal_stripping_radialness_ )
		return;
	_default_tidal_stripping_radialness_ = new_default_tidal_stripping_radialness;

	if(override_current)
			reset_tidal_stripping_radialness();
}
void brgastro::stripping_orbit::set_default_tidal_shocking_amplification(
		const double new_default_tidal_shocking_amplification)
{
	_default_tidal_shocking_amplification_ = new_default_tidal_shocking_amplification;
}
void brgastro::stripping_orbit::set_default_tidal_shocking_amplification(
		const double new_default_tidal_shocking_amplification,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_amplification == _default_tidal_shocking_amplification_ )
		return;

	if ( new_default_tidal_shocking_amplification < 0 )
	{
		throw std::logic_error("Attempt to set default tidal_shocking_amplification to value below minimum of 0.");
	}
	_default_tidal_shocking_amplification_ = new_default_tidal_shocking_amplification;

	if(override_current)
			reset_tidal_shocking_amplification();
}
void brgastro::stripping_orbit::set_default_tidal_shocking_persistance(
		const double new_default_tidal_shocking_persistance)
{
	_default_tidal_shocking_persistance_ = new_default_tidal_shocking_persistance;
}
void brgastro::stripping_orbit::set_default_tidal_shocking_persistance(
		const double new_default_tidal_shocking_persistance,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_persistance == _default_tidal_shocking_persistance_ )
		return;

	if ( new_default_tidal_shocking_persistance <= 0 )
	{
		throw std::logic_error("Attempt to set default tidal_stripping_persistance to value at or below minimum of 0.");
	}
	_default_tidal_shocking_persistance_ = new_default_tidal_shocking_persistance;

	if(override_current)
			reset_tidal_shocking_persistance();
}
void brgastro::stripping_orbit::set_default_tidal_shocking_power(
		const double new_default_tidal_shocking_power)
{
	_default_tidal_shocking_power_ = new_default_tidal_shocking_power;
}
void brgastro::stripping_orbit::set_default_tidal_shocking_power(
		const double new_default_tidal_shocking_power,
		const bool override_current,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_default_tidal_shocking_power == _default_tidal_shocking_power_ )
		return;
	_default_tidal_shocking_power_ = new_default_tidal_shocking_power;

	if(override_current)
			reset_tidal_shocking_power();
}
#endif

// Setting integration parameters
#if(1)
void brgastro::stripping_orbit::set_resolution( const size_t new_resolution,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_resolution == _base_resolution_ )
		return;

	if ( new_resolution < 2 )
	{
		throw std::logic_error("Attempt to set resolution to value below minimum of 2.");
	}
	clear_calcs();
	_base_resolution_ = new_resolution;
}
void brgastro::stripping_orbit::set_interpolation_type(
		const allowed_interpolation_type new_type,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( _interpolation_type_ == new_type )
		return;
	clear_calcs();
	_interpolation_type_ = new_type;
}
void brgastro::stripping_orbit::set_v_0( const double new_v_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_v_0 == _v_0_ )
		return;

	if ( new_v_0 <= 0 )
	{
		throw std::logic_error("Attempt to set v_0 to value at or below minimum of 0.");
	}
	clear_calcs();
	_v_0_ = new_v_0;
}
void brgastro::stripping_orbit::set_r_0( const double new_r_0,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_r_0 == _r_0_ )
		return;

	if ( new_r_0 <= 0 )
	{
		throw std::logic_error("Attempt to set r_0 to value at or below minimum of 0.");
	}
	clear_calcs();
	_r_0_ = new_r_0;
}
void brgastro::stripping_orbit::set_step_length_power( const double new_step_length_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_length_power == _step_length_power_ )
		return;

	clear_calcs();
	_step_length_power_ = new_step_length_power;
}
void brgastro::stripping_orbit::set_step_factor_max( const double new_step_factor_max,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_max == _step_factor_max_ )
		return;

	if ( new_step_factor_max < 1 )
	{
		throw std::logic_error("Attempt to set step_factor_max to value below minimum of 1.");
	}
	clear_calcs();
	_step_factor_max_ = new_step_factor_max;
}
void brgastro::stripping_orbit::set_step_factor_min( const double new_step_factor_min,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_factor_min == _step_factor_min_ )
		return;

	if ( new_step_factor_min > 1 )
	{
		throw std::logic_error("Attempt to set step_factor_min to value above minimum of 1.");
	}
	clear_calcs();
	_step_factor_min_ = new_step_factor_min;
}
#endif

// Setting tuning parameters
#if(1)

void brgastro::stripping_orbit::set_tidal_stripping_amplification(
		const double new_tidal_stripping_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_amplification == _tidal_stripping_amplification_ )
		return;

	if ( new_tidal_stripping_amplification < 0 )
	{
		throw std::logic_error("WARNING: Attempt to set tidal_stripping_amplification to value below minimum of 0.");
	}

	clear_calcs();
	_tidal_stripping_amplification_ = new_tidal_stripping_amplification;
}
void brgastro::stripping_orbit::set_tidal_stripping_deceleration(
		const double new_tidal_stripping_deceleration,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_deceleration == _tidal_stripping_deceleration_ )
		return;

	clear_calcs();
	_tidal_stripping_deceleration_ = new_tidal_stripping_deceleration;
}
void brgastro::stripping_orbit::set_tidal_stripping_radialness(
		const double new_tidal_stripping_radialness,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_radialness == _tidal_stripping_radialness_ )
		return;

	clear_calcs();
	_tidal_stripping_radialness_ = new_tidal_stripping_radialness;
}
void brgastro::stripping_orbit::set_tidal_shocking_amplification(
		const double new_tidal_shocking_amplification,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_amplification == _tidal_shocking_amplification_ )
		return;

	if ( new_tidal_shocking_amplification < 0 )
	{
		throw std::logic_error("WARNING: Attempt to set tidal_shocking_amplification to value below minimum of 0.");
	}

	clear_calcs();
	_tidal_shocking_amplification_ = new_tidal_shocking_amplification;
}
void brgastro::stripping_orbit::set_tidal_shocking_persistance(
		const double new_tidal_shocking_persistance,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_persistance == _tidal_shocking_persistance_ )
		return;

	if ( new_tidal_shocking_persistance <= 0 )
	{
		throw std::logic_error("WARNING: Attempt to set tidal_shocking_persistance to value at or below minimum of 0.");
	}
	clear_calcs();
	_tidal_shocking_persistance_ = new_tidal_shocking_persistance;
}
void brgastro::stripping_orbit::set_tidal_shocking_power(
		const double new_tidal_shocking_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_shocking_power == _tidal_shocking_power_ )
		return;

	clear_calcs();
	_tidal_shocking_power_ = new_tidal_shocking_power;
}
#endif

// Setting integration parameters
#if(1)
void brgastro::stripping_orbit::reset_resolution()
{
	// Check if anything is actually changing here
	if ( _base_resolution_ == _default_spline_resolution_ )
		return;

	clear_calcs();
	_base_resolution_ = _default_spline_resolution_;
}
void brgastro::stripping_orbit::reset_interpolation_type()
{
	// Check if anything is actually changing here
	if ( _interpolation_type_ == _default_interpolation_type_ )
		return;

	clear_calcs();
	_interpolation_type_ = _default_interpolation_type_;
}
void brgastro::stripping_orbit::reset_v_0()
{
	// Check if anything is actually changing here
	if ( _default_v_0_ == _v_0_ )
		return;

	clear_calcs();
	_v_0_ = _default_v_0_;
}
void brgastro::stripping_orbit::reset_r_0()
{
	// Check if anything is actually changing here
	if ( _default_r_0_ == _r_0_ )
		return;

	clear_calcs();
	_r_0_ = _default_r_0_;
}
void brgastro::stripping_orbit::reset_step_length_power()
{
	// Check if anything is actually changing here
	if ( _default_step_length_power_ == _step_length_power_ )
		return;

	clear_calcs();
	_step_length_power_ = _default_step_length_power_;
}
void brgastro::stripping_orbit::reset_step_factor_max()
{
	// Check if anything is actually changing here
	if ( _default_step_factor_max_ == _step_factor_max_ )
		return;

	clear_calcs();
	_step_factor_max_ = _default_step_factor_max_;
}
void brgastro::stripping_orbit::reset_step_factor_min()
{
	// Check if anything is actually changing here
	if ( _default_step_factor_min_ == _step_factor_min_ )
		return;

	clear_calcs();
	_step_factor_min_ = _default_step_factor_min_;
}
#endif

// Resetting tuning parameters
#if(1)

void brgastro::stripping_orbit::reset_tidal_stripping_amplification()
{
	// Check if anything is actually changing here
	if ( _default_tidal_stripping_amplification_ == _tidal_stripping_amplification_ )
		return;

	clear_calcs();
	_tidal_stripping_amplification_ = _default_tidal_stripping_amplification_;
}
void brgastro::stripping_orbit::reset_tidal_stripping_deceleration()
{
	// Check if anything is actually changing here
	if ( _default_tidal_stripping_deceleration_ == _tidal_stripping_deceleration_ )
		return;

	clear_calcs();
	_tidal_stripping_deceleration_ = _default_tidal_stripping_deceleration_;
}
void brgastro::stripping_orbit::reset_tidal_stripping_radialness()
{
	// Check if anything is actually changing here
	if ( _default_tidal_stripping_radialness_ == _tidal_stripping_radialness_ )
		return;

	clear_calcs();
	_tidal_stripping_radialness_ = _default_tidal_stripping_radialness_;
}
void brgastro::stripping_orbit::reset_tidal_shocking_amplification()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_amplification_ == _tidal_shocking_amplification_ )
		return;

	clear_calcs();
	_tidal_shocking_amplification_ = _default_tidal_shocking_amplification_;
}
void brgastro::stripping_orbit::reset_tidal_shocking_persistance()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_persistance_ == _tidal_shocking_persistance_ )
		return;

	clear_calcs();
	_tidal_shocking_persistance_ = _default_tidal_shocking_persistance_;
}
void brgastro::stripping_orbit::reset_tidal_shocking_power()
{
	// Check if anything is actually changing here
	if ( _default_tidal_shocking_power_ == _tidal_shocking_power_ )
		return;

	clear_calcs();
	_tidal_shocking_power_ = _default_tidal_shocking_power_;
}
#endif

void brgastro::stripping_orbit::set_t_min( CONST_BRG_TIME_REF new_t_min )
{
	_t_min_override_value_ = new_t_min;
	_override_t_min_ = true;
}

void brgastro::stripping_orbit::set_t_max( CONST_BRG_TIME_REF new_t_max )
{
	_t_max_override_value_ = new_t_max;
	_override_t_max_ = true;
}

void brgastro::stripping_orbit::reset_t_min()
{
	_t_min_natural_value_ = DBL_MAX;
	for ( size_t i = 0; i < _x_points_.size(); i++ )
	{
		if ( _x_points_.at( i ).first < _t_min_natural_value_ )
			_t_min_natural_value_ = _x_points_.at( i ).first;
	}
	_override_t_min_ = false;
}
void brgastro::stripping_orbit::reset_t_max()
{
	_t_max_natural_value_ = ( -DBL_MAX );
	for ( size_t i = 0; i < _x_points_.size(); i++ )
	{
		if ( _x_points_.at( i ).first > _t_max_natural_value_ )
			_t_max_natural_value_ = _x_points_.at( i ).first;
	}
	_override_t_max_ = false;
}

// Functions for determining how calc() will be called
void brgastro::stripping_orbit::set_record_full_data(
		const bool new_record_full_data ) const
{
	// Check if anything is actually changing here
	if ( new_record_full_data == _record_full_data_ )
		return;

	clear_calcs();
	_record_full_data_ = new_record_full_data;
}

// Function to calculate stripping
void brgastro::stripping_orbit::calc( const bool silent ) const
{
	std::vector< int > segment_resolutions( 0 );
	std::vector< bool > segments_to_skip( 0 );
	size_t num_segments_to_skip = 0;
	double t_min_to_use = (
			_override_t_min_ ? _t_min_override_value_ : _t_min_natural_value_ );
	double t_max_to_use = (
			_override_t_max_ ? _t_max_override_value_ : _t_max_natural_value_ );
	clear_calcs();

	// First, determine number of segments we'll be using

	// Start by generating the cleaned_discontinuity_times list
	_cleaned_discontinuity_times_.clear();
	for ( size_t i = 0; i < _discontinuity_times_.size(); i++ )
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

	for ( size_t i = 0; i < _x_points_.size(); i++ )
	{
		double t = _x_points_.at( i ).first;
		double x = _x_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_x_point( x, t );
	}

	for ( size_t i = 0; i < _y_points_.size(); i++ )
	{
		double t = _y_points_.at( i ).first;
		double y = _y_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_y_point( y, t );
	}

	for ( size_t i = 0; i < _z_points_.size(); i++ )
	{
		double t = _z_points_.at( i ).first;
		double z = _z_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_z_point( z, t );
	}

	for ( size_t i = 0; i < _vx_points_.size(); i++ )
	{
		double t = _vx_points_.at( i ).first;
		double vx = _vx_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vx_point( vx, t );
	}

	for ( size_t i = 0; i < _vy_points_.size(); i++ )
	{
		double t = _vy_points_.at( i ).first;
		double vy = _vy_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vy_point( vy, t );
	}

	for ( size_t i = 0; i < _vz_points_.size(); i++ )
	{
		double t = _vz_points_.at( i ).first;
		double vz = _vz_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_vz_point( vz, t );
	}

	for ( size_t i = 0; i < _vx_unknown_points_.size(); i++ )
	{
		double t = _vx_unknown_points_.at( i );
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vx_point( t );
	}

	for ( size_t i = 0; i < _vy_unknown_points_.size(); i++ )
	{
		double t = _vy_unknown_points_.at( i );
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vy_point( t );
	}

	for ( size_t i = 0; i < _vz_unknown_points_.size(); i++ )
	{
		double t = _vz_unknown_points_.at( i );
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_unknown_vz_point( t );
	}

	for ( size_t i = 0; i < _test_mass_points_.size(); i++ )
	{
		double t = _test_mass_points_.at( i ).first;
		double test_mass = _test_mass_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_test_mass_point( test_mass,
				t );
	}

	for ( size_t i = 0; i < _host_parameter_points_.size(); i++ )
	{
		double t = _host_parameter_points_.at( i ).first;
		std::vector< BRG_UNITS > host_parameters =
				_host_parameter_points_.at( i ).second;
		int segment_counter = 0;
		for ( size_t j = 0; j < _num_segments_ - 1; j++ )
		{
			if ( t > _cleaned_discontinuity_times_.at( j ) )
				segment_counter++;
		}
		_orbit_segments_.at( segment_counter ).add_host_parameter_point(
				host_parameters, t );
	}

	// Adjust t_min and t_max for each segment to leave no gaps, get segment resolution,
	// and check for segments with too few spline points
	if ( t_max_to_use <= t_min_to_use )
	{
		if(t_max_to_use==t_min_to_use)
			throw std::logic_error("t_min==t_max in stripping orbit; cannot integrate.");
		std::swap(t_max_to_use,t_min_to_use);
	}
	for ( size_t i = 0; i < _num_segments_; i++ )
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
		int segment_resolution = (int)( _base_resolution_
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
		throw std::runtime_error("Cannot calculate stripping for any segments of orbit.\n");
	}

	// Now loop through and calculate stripping for each segment in turn

	density_profile *temp_satellite = NULL;
	density_profile *temp_host = NULL;

	int last_good_segment = -1;

	try
	{
		for ( size_t i = 0; i < _num_segments_; i++ )
		{
			double frac_m_ret = 1;
			double fmvirret = 1;
			if ( !( segments_to_skip.at( i ) ) )
			{

				try {
					if ( last_good_segment == -1 ) // Special handling for first good segment
					{
						temp_satellite = _init_satellite_ptr_->density_profile_clone();
						temp_host = _init_host_ptr_->density_profile_clone();
					}
					else
					{
						delete temp_satellite;
						temp_satellite = NULL;
						delete temp_host;
						temp_host = NULL;
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
							throw std::runtime_error("ERROR: Could not connect orbit segments properly.\n");
						}
					}
					last_good_segment = i;

					// Pass parameters to the segment
					_pass_parameters_to_segment(_orbit_segments_.at( i ),
							temp_satellite,
							temp_host,
							segment_resolutions.at(i));
					// Calculate it at this stage, and check if it's disrupted or some other suspicious
					// result
					if (_orbit_segments_.at( i ).likely_disrupted())
					{
						_likely_disrupted_ = true;
						_bad_result_ = true;
						frac_m_ret = 0;
						fmvirret = 0;
						throw std::runtime_error("WARNING: Satellite halo likely disrupted.\n");
					}
					if (_orbit_segments_.at( i ).bad_result() )
					{
						_bad_result_ = true;
						throw std::runtime_error("WARNING: Could not calculate stripping for orbit segment.\n");
					}

					// If we're recording full data, fill up _m_ret_spline_points_ with new points
					if(_record_full_data_)
					{
						std::vector< std::pair< double, double > > new_m_rets =
								_orbit_segments_[i].m_ret_points();
						std::vector< std::pair< double, double > > new_m_vir_rets =
								_orbit_segments_[i].m_vir_ret_points();
						for(size_t j=0; j<new_m_rets.size(); j++)
						{
							try
							{
								_m_ret_interpolator_.try_add_point(
										new_m_rets[j].first,
										frac_m_ret*new_m_rets[j].second);
								_m_vir_ret_interpolator_.try_add_point(
										new_m_vir_rets[j].first,
										fmvirret*new_m_vir_rets[j].second);
							}
							catch(const std::exception &e)
							{
								std::cerr << "WARNING: Attempt to re-add point to m_ret_interpolator.\n";
								std::cerr.flush();
							}
						}
					}

					_orbit_segments_.at( i ).get_final_frac_m_ret( frac_m_ret );
					_orbit_segments_.at( i ).get_final_frac_m_vir_ret( fmvirret );

					// Record this as the final good segment, in case we don't find any more afterward
					_final_good_segment_ = _orbit_segments_.begin() + i;
				}
				catch (const std::exception &e) {
					if( !silent )
						std::cerr << e.what();
					if(_likely_disrupted_)
					{
						frac_m_ret = 0;
					}
					else
					{
						frac_m_ret = 1;
					}
					if ( _final_frac_m_ret_list_.size() > 0 )
					{
						_final_frac_m_ret_list_.push_back( frac_m_ret * _final_frac_m_ret_list_.back() );
						_final_frac_m_vir_ret_list_.push_back( fmvirret * _final_frac_m_vir_ret_list_.back() );
					}
					else
					{
						_final_frac_m_ret_list_.push_back( frac_m_ret );
						_final_frac_m_vir_ret_list_.push_back( fmvirret );
					}
					break;
				}

			}
			else
			{
				frac_m_ret = 1;
				fmvirret = 1;
			}

			if ( _final_frac_m_ret_list_.size() > 0 )
			{
				_final_frac_m_ret_list_.push_back( frac_m_ret * _final_frac_m_ret_list_.back() );
				_final_frac_m_vir_ret_list_.push_back( fmvirret * _final_frac_m_vir_ret_list_.back() );
			}
			else
			{
				_final_frac_m_ret_list_.push_back( frac_m_ret );
				_final_frac_m_vir_ret_list_.push_back( fmvirret );
			}
		}

		delete temp_satellite;
		temp_satellite = NULL;
		delete temp_host;
		temp_host = NULL;

		_calculated_ = true;
	}
	catch (const std::exception &e)
	{
		clear_calcs();
		delete temp_satellite;
		temp_satellite = NULL;
		delete temp_host;
		temp_host = NULL;

		_calculated_ = true;
		_bad_result_ = true;
	}
	return;

}

void brgastro::stripping_orbit::set_satellite_output_parameters(
		const std::vector< bool > & new_output_parameters )
{
	_satellite_output_parameters_ = new_output_parameters;
}
void brgastro::stripping_orbit::set_satellite_parameter_unitconvs(
		const std::vector< double > & new_parameter_unitconvs )
{
	_satellite_parameter_unitconvs_ = new_parameter_unitconvs;
}

void brgastro::stripping_orbit::set_host_output_parameters(
		const std::vector< bool > & new_output_parameters )
{
	_host_output_parameters_ = new_output_parameters;
}
void brgastro::stripping_orbit::set_host_parameter_unitconvs(
		const std::vector< double > & new_parameter_unitconvs )
{
	_host_parameter_unitconvs_ = new_parameter_unitconvs;
}

void brgastro::stripping_orbit::clear_satellite_output_parameters()
{
	_satellite_output_parameters_.clear();
}
void brgastro::stripping_orbit::clear_satellite_parameter_unitconvs()
{
	_satellite_parameter_unitconvs_.clear();
}

void brgastro::stripping_orbit::clear_host_output_parameters()
{
	_host_output_parameters_.clear();
}
void brgastro::stripping_orbit::clear_host_parameter_unitconvs()
{
	_host_parameter_unitconvs_.clear();
}

void brgastro::stripping_orbit::print_full_data( std::ostream *out ) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			std::cerr << "WARNING: Printing data for orbit which had a bad calculation.\n";
		}
	}

	// Loop through segments, printing each of them in turn
	bool first_good_segment = true;
	for ( size_t i = 0; i < _num_segments_; i++ )
	{
		if ( _orbit_segments_.at( i ).length() >= 2 )
		{
			_orbit_segments_.at( i ).set_satellite_output_parameters(
					_satellite_output_parameters_ );
			_orbit_segments_.at( i ).set_satellite_parameter_unitconvs(
					_satellite_parameter_unitconvs_ );
			_orbit_segments_.at( i ).set_host_output_parameters(
					_host_output_parameters_ );
			_orbit_segments_.at( i ).set_host_parameter_unitconvs(
					_host_parameter_unitconvs_ );
			if ( first_good_segment ) // Special handling for first segment
			{
				_orbit_segments_.at( i ).print_full_data( out, true, 1 );
				first_good_segment = false;
			}
			else
			{
				_orbit_segments_.at( i ).print_full_data( out, false,
						_final_frac_m_ret_list_.at( i - 1 ), _final_frac_m_vir_ret_list_.at( i - 1 ) );
				// Don't print header, m_ret_multiplier = frac_m_ret after last segment
			}
		}
	}
}

const int brgastro::stripping_orbit::get_final_m_ret( BRG_MASS & m_ret ) const
{
	if(_likely_disrupted_)
	{
		m_ret = 0;
		return 1;
	}
	else if(_bad_result_)
	{
		m_ret = _init_satellite_ptr_->mtot();
		return 1;
	}
	else if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			if(_likely_disrupted_)
			{
				m_ret = 0;
				return 1;
			}
			else
			{
				m_ret = _init_satellite_ptr_->mtot();
				return 1;
			}
		}
	}
	m_ret = _init_satellite_ptr_->mtot()*_final_frac_m_ret_list_.back();
	return 0;
}
const int brgastro::stripping_orbit::get_m_ret_at_t( CONST_BRG_TIME_REF t, BRG_MASS & m_ret) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		calc();
	}
	if( _bad_result_ )
	{
		m_ret = _init_satellite_ptr_->mtot();
		return 1;
	}
	m_ret = _init_satellite_ptr_->mtot()*max( _m_ret_interpolator_(t), 0. );
	return 0;
}
const int brgastro::stripping_orbit::get_final_frac_m_ret( double & frac_m_ret ) const
{
	if(likely_disrupted())
	{
		frac_m_ret = 0;
	}
	else
	{
		if ( !_calculated_ )
		{
			try
			{
				calc();
			}
			catch(const std::runtime_error &e)
			{
				return UNSPECIFIED_ERROR;
			}
		}
		if( _bad_result_ )
		{
			return UNSPECIFIED_ERROR;
		}
		frac_m_ret = _final_frac_m_ret_list_.back();
	}
	return 0;
}
const int brgastro::stripping_orbit::get_frac_m_ret_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	frac_m_ret = max( _m_ret_interpolator_(t), 0. );
	return 0;
}
const int brgastro::stripping_orbit::get_final_m_vir_ret( BRG_MASS & m_ret ) const
{
	if(likely_disrupted())
	{
		m_ret = 0;
	}
	else
	{
		if ( !_calculated_ )
		{
			try
			{
				calc();
			}
			catch(const std::runtime_error &e)
			{
				return UNSPECIFIED_ERROR;
			}
		}
		if( _bad_result_ )
		{
			return UNSPECIFIED_ERROR;
		}
		m_ret = _init_satellite_ptr_->mvir() *
				_final_frac_m_vir_ret_list_.back();
	}
	return 0;
}
const int brgastro::stripping_orbit::get_m_vir_ret_at_t( CONST_BRG_TIME_REF  t, BRG_MASS & m_ret) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	m_ret = _init_satellite_ptr_->mvir() *
			max( _m_vir_ret_interpolator_(t), 0. );
	return 0;
}
const int brgastro::stripping_orbit::get_final_frac_m_vir_ret( double & frac_m_ret ) const
{
	if(likely_disrupted())
	{
		frac_m_ret = 0;
	}
	else
	{
		if ( !_calculated_ )
		{
			try
			{
				calc();
			}
			catch(const std::runtime_error &e)
			{
				return UNSPECIFIED_ERROR;
			}
		}
		if( _bad_result_ )
		{
			return UNSPECIFIED_ERROR;
		}
		frac_m_ret = _final_frac_m_vir_ret_list_.back();
	}
	return 0;
}
const int brgastro::stripping_orbit::get_frac_m_vir_ret_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret) const
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	frac_m_ret = max( _m_vir_ret_interpolator_(t), 0. );
	return 0;
}
const int brgastro::stripping_orbit::get_comp_frac_m_ret_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret) const
{
	try
	{
		frac_m_ret = _test_mass_interpolator_(t);
	}
	catch(const std::exception &e)
	{
		return UNSPECIFIED_ERROR;
	}
	return 0;
}
const int brgastro::stripping_orbit::get_comp_frac_m_ret_error_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret) const
{
	try
	{
		frac_m_ret = _test_mass_error_interpolator_(t);
	}
	catch(const std::exception &e)
	{
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

const int brgastro::stripping_orbit::get_final_sum_deltarho(
		long double & final_sum_deltarho ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	_final_good_segment()->get_final_sum_deltarho( final_sum_deltarho );
	return 0;
}

const int brgastro::stripping_orbit::get_final_sum_deltarho(
		BRG_UNITS & final_sum_deltarho ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	_final_good_segment()->get_final_sum_deltarho( final_sum_deltarho );
	return 0;
}

const int brgastro::stripping_orbit::get_final_sum_gabdt(
		gabdt & final_sum_gabdt ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	_final_good_segment()->get_final_sum_gabdt( final_sum_gabdt );
	return 0;
}

const int brgastro::stripping_orbit::get_last_infall_time(
		BRG_TIME & t ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	t = _final_good_segment()->t_min_natural_value();
	return 0;
}

const int brgastro::stripping_orbit::clone_final_satellite(
		brgastro::density_profile * & final_satellite_clone ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
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
				final_satellite_clone = new brgastro::tNFW_profile;
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

const int brgastro::stripping_orbit::clone_final_host(
		brgastro::density_profile * & final_host_clone ) const
{
	if ( !_calculated_ )
	{
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
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
				final_host_clone = new brgastro::tNFW_profile;
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
const BRG_MASS brgastro::stripping_orbit::final_m_ret() const
{
	BRG_MASS result = -1;

	if(likely_disrupted()) return 0;

	if ( get_final_m_ret( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_m_ret.\n");
	}
	return result;
}
const BRG_MASS brgastro::stripping_orbit::m_ret_at_t(CONST_BRG_TIME_REF  t) const
{
	BRG_MASS result = -1;

	if ( get_m_ret_at_t( t, result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_m_ret.\n");
	}
	return result;
}
const BRG_MASS brgastro::stripping_orbit::final_m_vir_ret() const
{
	BRG_MASS result = -1;

	if(likely_disrupted()) return 0;

	if ( get_final_m_vir_ret( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_m_ret.\n");
	}
	return result;
}
const BRG_MASS brgastro::stripping_orbit::m_vir_ret_at_t(CONST_BRG_TIME_REF  t) const
{
	BRG_MASS result = -1;

	if ( get_m_vir_ret_at_t( t, result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_m_ret.\n");
	}
	return result;
}

const double brgastro::stripping_orbit::final_frac_m_vir_ret() const
{
	double result = -1;

	if ( get_final_frac_m_vir_ret( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_frac_m_ret.\n");
	}
	return result;
}
const double brgastro::stripping_orbit::frac_m_vir_ret_at_t(CONST_BRG_TIME_REF t) const
{
	double result = -1;

	if ( get_frac_m_vir_ret_at_t( t, result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_frac_m_ret.\n");
	}
	return result;
}
const double brgastro::stripping_orbit::comp_frac_m_ret_at_t(CONST_BRG_TIME_REF  t) const
{
	return _test_mass_interpolator_(t);
}
const double brgastro::stripping_orbit::comp_frac_m_ret_error_at_t(CONST_BRG_TIME_REF  t) const
{
	return _test_mass_error_interpolator_(t);
}

const BRG_UNITS brgastro::stripping_orbit::final_sum_deltarho() const
{
	BRG_UNITS result = -1;

	if ( get_final_sum_deltarho( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_sum_deltarho.\n");
	}
	return result;
}

const brgastro::gabdt brgastro::stripping_orbit::final_sum_gabdt() const
{
	brgastro::gabdt result;

	if ( get_final_sum_gabdt( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::final_sum_gabdt.\n");
	}
	return result;
}

const BRG_TIME brgastro::stripping_orbit::last_infall_time() const
{
	BRG_TIME result;

	if ( get_last_infall_time( result ) )
	{
		throw std::runtime_error("ERROR: Could not calculate stripping in stripping_orbit::last_infall_time.\n");
	}
	return result;
}

const bool brgastro::stripping_orbit::likely_disrupted() const
{
	if(!_calculated_)
	{
		try
		{
			calc();
		}
		catch(const std::exception &e)
		{
			// Do nothing on exception here
		}
	}
	return _likely_disrupted_;
}

const brgastro::density_profile * brgastro::stripping_orbit::final_satellite() const
{
	if ( !_calculated_ )
	{
		if(_satellite_loaded_)
		{
			try
			{
				calc();
			}
			catch(const std::runtime_error &e)
			{
				return _init_satellite_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_satellite() without init_satellite assigned.\n");
			return NULL;
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) || (_likely_disrupted_) )
		return _init_satellite_ptr_; // Sanest/safest option in this case
	else
		return _final_good_segment()->final_satellite();
}

const brgastro::density_profile * brgastro::stripping_orbit::final_host() const
{
	if ( !_calculated_ )
	{
		if(_host_loaded_)
		{
			try
			{
				calc();
			}
			catch(const std::runtime_error &e)
			{
				return _init_host_ptr_;
			}
		}
		else
		{
			throw std::runtime_error("ERROR: Attempt to call stripping_orbit::final_host() without init_host assigned.\n");
			return NULL;
		}
	}
	if( (_bad_result_) || (_final_good_segment()==_orbit_segments_.end()) || (_likely_disrupted_) )
		return _init_host_ptr_; // Sanest/safest option in this case
	else
		return _final_good_segment()->final_host();
}

const int brgastro::stripping_orbit::get_quality_of_fit( double & Q, const bool use_virial,
		const unsigned int samples)
{
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		set_record_full_data( true );
		try
		{
			calc();
		}
		catch(const std::runtime_error &e)
		{
			return UNSPECIFIED_ERROR;
		}
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}
	double temp_Q = 0;
	double t_step = (t_max()-t_min())/samples;

	interpolator * interp=NULL;
	if(use_virial)
		interp = &_m_vir_ret_interpolator_;
	else
		interp = &_m_ret_interpolator_;



	for(double t=t_min(), tm=t_max(); t<tm; t+=t_step)
	{
		temp_Q += square(((*interp)(t) - _test_mass_interpolator_(t))
				/safe_d( _test_mass_error_interpolator_(t) ));
	}

	Q = temp_Q/safe_d(samples);

	return 0;
}

const double brgastro::stripping_orbit::quality_of_fit(const bool use_virial, const unsigned int samples)
{
	double Q=-1;
	if(get_quality_of_fit(Q,use_virial,samples))
	{
		throw std::runtime_error("Cannot determine quality of fit for stripping_orbit.");
	}
	return Q;
}

#endif // end brgastro::stripping_orbit_segment class function definitions


#endif // end class function definitions
