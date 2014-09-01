/**********************************************************************\
  @file stripping_orbit_segment.cpp

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

#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "brg/global.h"

#include "brg/file_functions.h"

#include "brg/math/interpolator/interpolator.h"
#include "brg/math/interpolator/interpolator_derivative.h"

#include "brg/physics/astro/SALTSA/gabdt.h"

#include "stripping_orbit_segment.h"


// brgastro::stripping_orbit_segment class method implementations
#if (1)

brgastro::stripping_orbit_segment::stripping_orbit_segment()
{
	_init();
}

const int brgastro::stripping_orbit_segment::_init()
{
	_current_satellite_in_use_ = false;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	return clear();
}

// Swap functions
void brgastro::stripping_orbit_segment::swap(stripping_orbit_segment &other)
{
	using std::swap;

	// Integration parameters
#if(1)
	swap(_spline_resolution_, other._spline_resolution_);
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

	swap(_rt_list_, other._rt_list_);
	swap(_rt_ratio_list_, other._rt_ratio_list_);
	swap(_delta_rho_list_, other._delta_rho_list_);
	swap(_sum_delta_rho_list_, other._sum_delta_rho_list_);
	swap(_x_data_, other._x_data_);
	swap(_y_data_, other._y_data_);
	swap(_z_data_, other._z_data_);
	swap(_vx_data_, other._vx_data_);
	swap(_vy_data_, other._vy_data_);
	swap(_vz_data_, other._vz_data_);
	swap(_t_min_natural_value_, other._t_min_natural_value_);
	swap(_t_max_natural_value_, other._t_max_natural_value_);
	swap(_satellite_parameter_data_, other._satellite_parameter_data_);
	swap(_host_parameter_data_, other._host_parameter_data_);
	swap(_num_parameters_, other._num_parameters_);
	swap(_m_ret_list_, other._m_ret_list_);
	swap(_m_vir_ret_list_, other._m_vir_ret_list_);
	swap(_satellite_parameter_unitconvs_,
			other._satellite_parameter_unitconvs_);
	swap(_satellite_output_parameters_,
			other._satellite_output_parameters_);
	swap(_host_parameter_unitconvs_, other._host_parameter_unitconvs_);
	swap(_host_output_parameters_, other._host_output_parameters_);
	swap(_gabdt_list_, other._gabdt_list_);
	swap(_sum_gabdt_list_, other._sum_gabdt_list_);
	swap(_phase_list_, other._phase_list_);
	swap(_phase_output_list_, other._phase_output_list_);
	swap(_record_full_data_, other._record_full_data_);
	swap(_x_spline_, other._x_spline_);
	swap(_y_spline_, other._y_spline_);
	swap(_z_spline_, other._z_spline_);
	swap(_vx_spline_, other._vx_spline_);
	swap(_vy_spline_, other._vy_spline_);
	swap(_vz_spline_, other._vz_spline_);
	swap(_test_mass_spline_, other._test_mass_spline_);
	swap(_host_parameter_splines_, other._host_parameter_splines_);
	swap(_init_host_ptr_, other._init_host_ptr_);
	swap(_init_satellite_ptr_, other._init_satellite_ptr_);
	swap(_init_sum_delta_rho_, other._init_sum_delta_rho_);
	swap(_init_sum_gabdt_, other._init_sum_gabdt_);
	swap(_calculated_, other._calculated_);
	swap(_bad_result_, other._bad_result_);
	swap(_current_satellite_in_use_, other._current_satellite_in_use_);
	swap(_current_satellite_ptr_, other._current_satellite_ptr_);
	swap(_current_host_in_use_, other._current_host_in_use_);
	swap(_current_host_ptr_, other._current_host_ptr_);
	swap(_host_loaded_, other._host_loaded_);
	swap(_satellite_loaded_, other._satellite_loaded_);
	swap(_using_private_init_host_, other._using_private_init_host_);
	swap(_using_private_init_satellite_,
			other._using_private_init_satellite_);
	swap(_private_tNFW_init_host_, other._private_tNFW_init_host_);
	swap(_private_tNFW_init_satellite_,
			other._private_tNFW_init_satellite_);
	swap(_evolving_host_, other._evolving_host_);
	swap(_override_t_max_, other._override_t_max_);
	swap(_override_t_min_, other._override_t_min_);
	swap(_t_max_override_val_, other._t_max_override_val_);
	swap(_t_min_override_val_, other._t_min_override_val_);
	swap(_likely_disrupted_, other._likely_disrupted_);

	if ( _using_private_init_host_ )
	{
		_init_host_ptr_ = &_private_tNFW_init_host_;
	}
	if ( other._using_private_init_host_ )
	{
		other._init_host_ptr_ = &other._private_tNFW_init_host_;
	}
	if ( _using_private_init_satellite_ )
	{
		_init_satellite_ptr_ = &_private_tNFW_init_satellite_;
	}
	if ( other._using_private_init_satellite_ )
	{
		other._init_satellite_ptr_ = &other._private_tNFW_init_satellite_;
	}
}
namespace std
{
	template <>
	void swap(brgastro::stripping_orbit_segment &same,
			brgastro::stripping_orbit_segment &other)
	{
		same.swap(other);
	}
}

brgastro::stripping_orbit_segment::stripping_orbit_segment(
		const stripping_orbit_segment &other )
{

	// Integration parameters
#if(1)
	_spline_resolution_ = other._spline_resolution_;
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

	_rt_list_ = other._rt_list_;
	_rt_ratio_list_ = other._rt_ratio_list_;
	_delta_rho_list_ = other._delta_rho_list_;
	_sum_delta_rho_list_ = other._sum_delta_rho_list_;
	_x_data_ = other._x_data_;
	_y_data_ = other._y_data_;
	_z_data_ = other._z_data_;
	_vx_data_ = other._vx_data_;
	_vy_data_ = other._vy_data_;
	_vz_data_ = other._vz_data_;
	_t_min_natural_value_ = other._t_min_natural_value_;
	_t_max_natural_value_ = other._t_max_natural_value_;
	_satellite_parameter_data_ = other._satellite_parameter_data_;
	_host_parameter_data_ = other._host_parameter_data_;
	_num_parameters_ = other._num_parameters_;
	_m_vir_ret_list_ = other._m_vir_ret_list_;
	_satellite_parameter_unitconvs_ =
			other._satellite_parameter_unitconvs_;
	_satellite_output_parameters_ =
			other._satellite_output_parameters_;
	_host_parameter_unitconvs_ = other._host_parameter_unitconvs_;
	_host_output_parameters_ = other._host_output_parameters_;
	_gabdt_list_ = other._gabdt_list_;
	_sum_gabdt_list_ = other._sum_gabdt_list_;
	_phase_list_ = other._phase_list_;
	_phase_output_list_ = other._phase_output_list_;
	_record_full_data_ = other._record_full_data_;
	_x_spline_ = other._x_spline_;
	_y_spline_ = other._y_spline_;
	_z_spline_ = other._z_spline_;
	_vx_spline_ = other._vx_spline_;
	_vy_spline_ = other._vy_spline_;
	_vz_spline_ = other._vz_spline_;
	_test_mass_spline_ = other._test_mass_spline_;
	_host_parameter_splines_ = other._host_parameter_splines_;
	_init_host_ptr_ = other._init_host_ptr_;
	_init_satellite_ptr_ = other._init_satellite_ptr_;
	_init_sum_delta_rho_ = other._init_sum_delta_rho_;
	_init_sum_gabdt_ = other._init_sum_gabdt_;
	_calculated_ = other._calculated_;
	_bad_result_ = other._bad_result_;
	_current_satellite_in_use_ = other._current_satellite_in_use_;
	_current_host_in_use_ = other._current_host_in_use_;
	_host_loaded_ = other._host_loaded_;
	_satellite_loaded_ = other._satellite_loaded_;
	_using_private_init_host_ = other._using_private_init_host_;
	_using_private_init_satellite_ =
			other._using_private_init_satellite_;
	_private_tNFW_init_host_ = other._private_tNFW_init_host_;
	_private_tNFW_init_satellite_ =
			other._private_tNFW_init_satellite_;
	_evolving_host_ = other._evolving_host_;
	_override_t_max_ = other._override_t_max_;
	_override_t_min_ = other._override_t_min_;
	_t_max_override_val_ = other._t_max_override_val_;
	_t_min_override_val_ = other._t_min_override_val_;
	_likely_disrupted_ = other._likely_disrupted_;

	if ( _current_satellite_in_use_ )
	{
		_current_satellite_ptr_ =
				other._current_satellite_ptr_->density_profile_clone();
	}
	else
	{
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		_current_host_ptr_ =
				other._current_host_ptr_->density_profile_clone();
	}
	else
	{
		_current_host_ptr_ = NULL;
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

brgastro::stripping_orbit_segment & brgastro::stripping_orbit_segment::operator=(
		stripping_orbit_segment other )
{
	swap(other);
	return *this;
}

brgastro::stripping_orbit_segment::stripping_orbit_segment(
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

brgastro::stripping_orbit_segment::~stripping_orbit_segment()
{
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
}

brgastro::stripping_orbit_segment *brgastro::stripping_orbit_segment::stripping_orbit_spline_clone() const
{
	return new stripping_orbit_segment( *this );
}

const int brgastro::stripping_orbit_segment::clear()
{
	clear_calcs();

	// Integration parameters
#if(1)
	_spline_resolution_ = brgastro::stripping_orbit::default_spline_resolution();
	_interpolation_type_ = brgastro::stripping_orbit::default_interpolation_type();
	_v_0_ = brgastro::stripping_orbit::default_spline_resolution();
	_r_0_ = brgastro::stripping_orbit::default_spline_resolution();
	_step_length_power_ = brgastro::stripping_orbit::default_spline_resolution();
	_step_factor_max_ = brgastro::stripping_orbit::default_spline_resolution();
	_step_factor_min_ = brgastro::stripping_orbit::default_spline_resolution();
#endif

	// Tuning values
#if(1)
	_tidal_stripping_amplification_ = brgastro::stripping_orbit::default_tidal_stripping_amplification();
	_tidal_stripping_deceleration_ = brgastro::stripping_orbit::default_tidal_stripping_deceleration();
	_tidal_stripping_radialness_ = brgastro::stripping_orbit::default_tidal_stripping_radialness();
	_tidal_shocking_amplification_ = brgastro::stripping_orbit::default_tidal_shocking_amplification();
	_tidal_shocking_persistance_ = brgastro::stripping_orbit::default_tidal_shocking_persistance();
	_tidal_shocking_power_ = brgastro::stripping_orbit::default_tidal_shocking_power();
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

	_init_host_ptr_ = NULL;
	_init_satellite_ptr_ = NULL;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}

	_current_satellite_ptr_ = NULL;
	_current_satellite_in_use_ = false;
	_current_host_ptr_ = NULL;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_satellite_loaded_ = false;
	_using_private_init_host_ = false;
	_using_private_init_satellite_ = false;
	_evolving_host_ = false;

	return 0;
}

const int brgastro::stripping_orbit_segment::clear_calcs() const
{
	_rt_list_.clear();
	_rt_ratio_list_.clear();
	_phase_output_list_.clear();
	_m_ret_list_.clear();
	_m_vir_ret_list_.clear();
	_delta_rho_list_.clear();
	_sum_delta_rho_list_.clear();
	_gabdt_list_.clear();
	_sum_gabdt_list_.clear();

	_calculated_ = false;
	_bad_result_ = false;
	_likely_disrupted_ = false;
	return 0;
}

// Setting integration parameters
#if(1)
const int brgastro::stripping_orbit_segment::set_resolution( const int new_resolution,
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
const int brgastro::stripping_orbit_segment::set_interpolation_type(
		const brgastro::stripping_orbit::allowed_interpolation_type new_type,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( _interpolation_type_ == new_type )
		return 0;
	clear_calcs();
	_interpolation_type_ = new_type;
	return 0;
}
const int brgastro::stripping_orbit_segment::set_v_0( const double new_v_0,
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
const int brgastro::stripping_orbit_segment::set_r_0( const double new_r_0,
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
const int brgastro::stripping_orbit_segment::set_step_length_power( const double new_step_length_power,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_step_length_power == _step_length_power_ )
		return 0;

	clear_calcs();
	_step_length_power_ = new_step_length_power;
	return 0;
}
const int brgastro::stripping_orbit_segment::set_step_factor_max( const double new_step_factor_max,
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
const int brgastro::stripping_orbit_segment::set_step_factor_min( const double new_step_factor_min,
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

const int brgastro::stripping_orbit_segment::set_tidal_stripping_amplification(
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
const int brgastro::stripping_orbit_segment::set_tidal_stripping_deceleration(
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
const int brgastro::stripping_orbit_segment::set_tidal_stripping_radialness(
		const double new_tidal_stripping_radialness,
		const bool silent )
{
	// Check if anything is actually changing here
	if ( new_tidal_stripping_radialness == _tidal_stripping_radialness_ )
		return 0;

	clear_calcs();
	_tidal_stripping_radialness_ = new_tidal_stripping_radialness;
	return 0;
}
const int brgastro::stripping_orbit_segment::set_tidal_shocking_amplification(
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
const int brgastro::stripping_orbit_segment::set_tidal_shocking_persistance(
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
const int brgastro::stripping_orbit_segment::set_tidal_shocking_power(
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


const int brgastro::stripping_orbit_segment::add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
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

const int brgastro::stripping_orbit_segment::add_point( CONST_BRG_DISTANCE_REF x,
		CONST_BRG_DISTANCE_REF y, CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
		CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
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

const int brgastro::stripping_orbit_segment::add_x_point(
		CONST_BRG_DISTANCE_REF x, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_y_point(
		CONST_BRG_DISTANCE_REF y, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_z_point(
		CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_vx_point(
		CONST_BRG_VELOCITY_REF vx, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_vy_point(
		CONST_BRG_VELOCITY_REF vy, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_vz_point(
		CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_unknown_vx_point(
		CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_unknown_vy_point(
		CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_unknown_vz_point(
		CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_test_mass_point(
		const double test_mass, CONST_BRG_TIME_REF t )
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

const int brgastro::stripping_orbit_segment::add_host_parameter_point(
		const unsigned int num_parameters,
		const std::vector< BRG_UNITS > & parameters, CONST_BRG_TIME_REF t,
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

const int brgastro::stripping_orbit_segment::_reserve( const int n,
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
	_m_ret_list_.reserve( n );
	_m_vir_ret_list_.reserve( n );
	_phase_output_list_.reserve( n );
	_rt_list_.reserve( n );
	_rt_ratio_list_.reserve( n );
	_delta_rho_list_.reserve( n );
	_sum_delta_rho_list_.reserve( n );
	_gabdt_list_.reserve( n );
	_sum_gabdt_list_.reserve( n );
	return 0;
}

const int brgastro::stripping_orbit_segment::set_init_host(
		const density_profile *new_init_host )
{
	_init_host_ptr_ = new_init_host;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
	_current_host_ptr_ = _init_host_ptr_->density_profile_clone();
	_current_host_in_use_ = true;
	_host_loaded_ = true;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_init_satellite(
		const density_profile *new_init_satellite )
{
	_init_satellite_ptr_ = new_init_satellite;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	_current_satellite_ptr_ = _init_satellite_ptr_->density_profile_clone();
	_current_satellite_in_use_ = true;
	_satellite_loaded_ = true;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_init_sum_deltarho(
		const long double &new_init_sum_deltarho )
{
	_init_sum_delta_rho_ = new_init_sum_deltarho;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_init_sum_gabdt(
		const gabdt &new_init_sum_gabdt )
{
	_init_sum_gabdt_ = new_init_sum_gabdt;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_tNFW_init_satellite(
		CONST_BRG_MASS_REF new_init_mvir0, const double z,
		const double new_init_c, const double new_init_tau )
{
	_using_private_init_satellite_ = true;
	_private_tNFW_init_satellite_ = tNFW_profile( new_init_mvir0, z,
			new_init_c, new_init_tau );
	_init_satellite_ptr_ = &_private_tNFW_init_satellite_;

	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
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

const int brgastro::stripping_orbit_segment::set_tNFW_host(
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

			if ( _current_satellite_in_use_ )
			{
				delete _current_satellite_ptr_;
				_current_satellite_ptr_ = NULL;
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

const int brgastro::stripping_orbit_segment::set_t_min(
		CONST_BRG_TIME_REF new_t_min )
{
	_t_min_override_val_ = new_t_min;
	_override_t_min_ = true;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_t_max(
		CONST_BRG_TIME_REF new_t_max )
{
	_t_max_override_val_ = new_t_max;
	_override_t_max_ = true;
	return 0;
}

const int brgastro::stripping_orbit_segment::reset_t_min()
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

const int brgastro::stripping_orbit_segment::reset_t_max()
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

const int brgastro::stripping_orbit_segment::clear_points()
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
const int brgastro::stripping_orbit_segment::clear_host_parameter_points()
{
	_host_parameter_splines_.clear();
	_calculated_ = false;
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_init_sum_deltarho()
{
	_init_sum_delta_rho_ = 0;
	_calculated_ = false;
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_init_sum_gabdt()
{
	_init_sum_gabdt_.clear();
	_init_sum_gabdt_.override_zero();
	_calculated_ = false;
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_init_satellite()
{
	_init_satellite_ptr_ = NULL;
	if ( _current_satellite_in_use_ )
	{
		delete _current_satellite_ptr_;
		_current_satellite_ptr_ = NULL;
	}
	_current_satellite_ptr_ = 0;
	_current_satellite_in_use_ = false;
	_satellite_loaded_ = false;
	_using_private_init_satellite_ = false;
	_calculated_ = false;
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_init_host()
{
	_init_host_ptr_ = NULL;
	if ( _current_host_in_use_ )
	{
		delete _current_host_ptr_;
		_current_host_ptr_ = NULL;
	}
	_current_host_ptr_ = NULL;
	_current_host_in_use_ = false;
	_host_loaded_ = false;
	_using_private_init_host_ = false;
	_calculated_ = false;
	return 0;
}

const unsigned int brgastro::stripping_orbit_segment::length() const
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
const int brgastro::stripping_orbit_segment::set_record_full_data(
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
const int brgastro::stripping_orbit_segment::calc( const bool silent ) const
{
	BRG_DISTANCE r;
	BRG_VELOCITY v, vt, vr;
	BRG_TIME t, t_step, t_min_to_use, t_max_to_use;
	BRG_UNITS current_deltarho;
	std::vector< BRG_UNITS > satellite_parameters;
	std::vector< BRG_UNITS > host_parameters;
	int num_host_parameters = 0;
	int counter = 0;
	double step_length_factor = 1;
	double m_ret = 1;

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

	// Set up interpolation type
	_pass_interpolation_type();

	_m_ret_list_.push_back( 1 );
	_m_vir_ret_list_.push_back( 1 );
	_delta_rho_list_.push_back( current_deltarho );
	_sum_delta_rho_list_.push_back( _init_sum_delta_rho_ );

	_sum_gabdt_list_.push_back( current_gabdt );
	current_gabdt.set_host_ptr( _init_host_ptr_ );
	current_gabdt.override_zero();
	_gabdt_list_.push_back( current_gabdt );
	current_gabdt.override_zero();

	// Calculate radius and velocity at init position
	r = brgastro::dist3d( _x_spline_( t_min_to_use ),
			_y_spline_( t_min_to_use ), _z_spline_( t_min_to_use ) );
	v = brgastro::dist3d( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ) );
	vr = brgastro::dot_product( _vx_spline_( t_min_to_use ),
			_vy_spline_( t_min_to_use ), _vz_spline_( t_min_to_use ),
			_x_spline_( t_min_to_use ), _y_spline_( t_min_to_use ),
			_z_spline_( t_min_to_use ) ) / r;
	if ( std::fabs( vr ) >= std::fabs( v ) )
		vt = 0;
	else
		vt = brgastro::quad_sub( v, vr );

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
			_get_rt( _current_host_ptr_, _current_satellite_ptr_, r,
					vr, vt, t_step, current_deltarho ) );
	_rt_ratio_list_.push_back( _rt_list_.at( 0 ) / _rvir() );

	// Get satellite parameters at init position
	if ( _record_full_data_ )
	{
		satellite_parameters = _current_satellite_ptr_->get_parameters( );
		_satellite_parameter_data_.push_back( satellite_parameters );
		host_parameters = _current_host_ptr_->get_parameters( );
		_host_parameter_data_.push_back( host_parameters );
	}

	const BRG_MASS init_mtot = _current_satellite_ptr_->mtot();
	const BRG_MASS init_mvir = _current_satellite_ptr_->mvir();

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
			r = brgastro::dist3d( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			v = brgastro::quad_add( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ) );
			vr = brgastro::dot_product( _vx_spline_( t ), _vy_spline_( t ),
					_vz_spline_( t ), _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) ) / r;
			if ( std::fabs( vr ) >= std::fabs( v ) )
				vt = 0;
			else
				vt = brgastro::quad_sub( v, vr );
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
				_current_host_ptr_->set_parameters( host_parameters );
			}

			// Calculate effects of tidal stripping and shocks
			// Effects of stripping
			m_ret = _tidal_strip_retained( _current_host_ptr_,
					_current_satellite_ptr_, r, vr, vt,
					t_step * step_length_factor,
					_sum_delta_rho_list_.at( counter - 1 ) );
			_m_ret_list_.push_back( _m_ret_list_.at( counter - 1 ) * m_ret );

			// Calculate adjusted fraction (so numerical errors don't cause a consistent
			// offset to add up).

			double adjusted_m_ret = _m_ret_list_.back() * init_mtot/safe_d(_current_satellite_ptr_->mtot());
			if( isbad(adjusted_m_ret) or (adjusted_m_ret>1.1) or (adjusted_m_ret<0)) adjusted_m_ret = 0;

			_current_satellite_ptr_->truncate_to_fraction( adjusted_m_ret );

			_m_vir_ret_list_.push_back( _current_satellite_ptr_->mvir() /
					init_mvir);

			if(m_ret == 0)
			{
				_likely_disrupted_ = true;
				_bad_result_ = true;
				break;
			}

			// Effects of shocking

			double t_shock = r / safe_d(v);
			double t_recover = _current_satellite_ptr_->othmtot()/(2*pi);
			double x;
			double gabdt_scaling_factor;
			x = max(0.,t_shock / safe_d(t_recover) / ( 2 * pi ));

			if(isbad(t_recover) or (t_recover<=0) or isbad(x) or (x<=0))
			{
				_likely_disrupted_ = true;
				_bad_result_ = true;
				break;
			}
			else
			{
				gabdt_scaling_factor = max(0.,1
									- ( t_step * step_length_factor
											/ safe_d( t_recover * _tidal_shocking_persistance_ ) ) );
			}

			current_gabdt.set_host_ptr( _current_host_ptr_ );
			current_gabdt.set_pos( _x_spline_( t ), _y_spline_( t ),
					_z_spline_( t ) );
			current_gabdt.set_dt( t_step * step_length_factor );
			if(current_gabdt.calc_dv())
			{
				_bad_result_ = true;
				break;
			}
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
					_get_rt( _current_host_ptr_,
							_current_satellite_ptr_, r, vr, vt,
							t_step * step_length_factor,
							_sum_delta_rho_list_.at( counter - 1 ) ) );
			_rt_ratio_list_.push_back( _rt_list_.at( counter ) / _rvir( 0 ) );

			satellite_parameters = _current_satellite_ptr_->get_parameters(  );
			_satellite_parameter_data_.push_back( satellite_parameters );

			host_parameters =_current_host_ptr_->get_parameters(  );
			_host_parameter_data_.push_back( host_parameters );

		}
		catch ( const std::exception & e )
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

	if((_bad_result_) || (_likely_disrupted_)) return UNSPECIFIED_ERROR;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_satellite_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_satellite_output_parameters )
{
	_satellite_output_parameters_ = new_satellite_output_parameters;
	return 0;
}
const int brgastro::stripping_orbit_segment::set_satellite_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > &new_satellite_parameter_unitconvs )
{
	_satellite_parameter_unitconvs_ = new_satellite_parameter_unitconvs;
	return 0;
}

const int brgastro::stripping_orbit_segment::set_host_output_parameters(
		const unsigned int num_out_parameters,
		const std::vector< bool > & new_host_output_parameters )
{
	_host_output_parameters_ = new_host_output_parameters;
	return 0;
}
const int brgastro::stripping_orbit_segment::set_host_parameter_unitconvs(
		const unsigned int num_out_parameters,
		const std::vector< double > & new_host_parameter_unitconvs )
{
	_host_parameter_unitconvs_ = new_host_parameter_unitconvs;
	return 0;
}

const int brgastro::stripping_orbit_segment::clear_satellite_output_parameters()
{
	_satellite_output_parameters_.clear();
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_satellite_parameter_unitconvs()
{
	_satellite_parameter_unitconvs_.clear();
	return 0;
}

const int brgastro::stripping_orbit_segment::clear_host_output_parameters()
{
	_host_output_parameters_.clear();
	return 0;
}
const int brgastro::stripping_orbit_segment::clear_host_parameter_unitconvs()
{
	_host_parameter_unitconvs_.clear();
	return 0;
}

const int brgastro::stripping_orbit_segment::print_full_data(
		std::ostream *out, const bool include_header,
		const double m_ret_multiplier, const double m_vir_ret_multiplier, const bool silent ) const
{
	const int num_columns_base = 18;
	int num_extra_satellite_columns = 0;
	int num_extra_host_columns = 0;
	int num_rows = 0;
	std::vector< std::string > header;
	std::vector< std::vector< std::string > > data;
	std::stringstream ss;

	// Calculate data if necessary
	if ( ( !_calculated_ ) || ( !_record_full_data_ ) )
	{
		if ( int errcode = set_record_full_data( true ) )
			return errcode + LOWER_LEVEL_ERROR;
		if ( calc() )
			return UNSPECIFIED_ERROR;
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
			catch ( const std::exception &e )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading satellite_output_parameters.\n"
							<< "No extra parameters will be output.\n"
							<< e.what();
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
			catch ( const std::exception &e )
			{
				if ( !silent )
					std::cerr
							<< "WARNING: Error reading host_output_parameters.\n"
							<< "No extra parameters will be output.\n"
							<< e.what();
				num_extra_host_columns = 0;
			}
		}
	}

	int num_columns = num_columns_base + num_extra_satellite_columns
			+ num_extra_host_columns;

	header.resize( num_columns );
	make_array2d( data, num_columns, num_rows );

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
	header[12] = "m_vir_ret";
	header[13] = "m_vir_frac_lost";
	header[14] = "comp_m_ret";
	header[15] = "comp_m_frac_lost";
	header[16] = "rt";
	header[17] = "rt/rvir";

	if ( num_extra_satellite_columns > 0 )
	{
		int extra_column_counter = 0;
		std::vector< std::string > parameter_names( 0 );

		parameter_names = _current_satellite_ptr_->get_parameter_names(  );

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
				catch ( const std::out_of_range & )
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
		std::vector< std::string > parameter_names( 0 );

		parameter_names = _current_host_ptr_->get_parameter_names(  );

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
				catch ( const std::out_of_range & )
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

	BRG_TIME t_min_to_use, t_max_to_use;

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
				<< brgastro::dist3d( _phase_output_list_.at( i ).x,
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
		ss << _m_ret_list_.at( i ) * m_ret_multiplier;
		data[10][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss << ( 1 - _m_ret_list_.at( i ) / _m_ret_list_.at( i - 1 ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << ( 1 - _m_ret_list_.at( i+1 ) / _m_ret_list_.at( i ) )
							/ ( _phase_output_list_.at( i+1 ).t
									- _phase_output_list_.at( i ).t );
		data[11][i] = ss.str();
		ss.str( "" );
		ss << _m_vir_ret_list_.at( i ) * m_vir_ret_multiplier;
		data[12][i] = ss.str();
		ss.str( "" );
		if ( i > 0 )
			ss << ( 1 - _m_vir_ret_list_.at( i ) / _m_vir_ret_list_.at( i - 1 ) )
							/ ( _phase_output_list_.at( i ).t
									- _phase_output_list_.at( i - 1 ).t );
		else
			ss << ( 1 - _m_vir_ret_list_.at( i+1 ) / _m_vir_ret_list_.at( i ) )
							/ ( _phase_output_list_.at( i+1 ).t
									- _phase_output_list_.at( i ).t );
		data[13][i] = ss.str();
		ss.str( "" );
		ss << _test_mass_spline_( _phase_output_list_.at( i ).t );
		data[14][i] = ss.str();
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
		data[15][i] = ss.str();
		ss.str( "" );
		ss << _rt_list_.at( i ) * unitconv::mtokpc;
		data[16][i] = ss.str();
		ss.str( "" );
		ss << _rt_ratio_list_.at( i );
		data[17][i] = ss.str();

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
					catch ( const std::out_of_range & )
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
					catch ( const std::out_of_range & )
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
		print_table( *out, data, header, silent );
	else
		print_table( *out, data, std::vector<std::string>(), silent );
	return 0;
}

const BRG_UNITS brgastro::stripping_orbit_segment::_delta_rho(
		const int index, const double x, CONST_BRG_TIME_REF t_step,
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
			std::cerr << "ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n"
					<< "x = " << x << std::endl
					<< "_tidal_shocking_power_ = " << _tidal_shocking_power_ << std::endl
					<< "_tidal_shocking_amplification_ = " << _tidal_shocking_amplification_ << std::endl
					<< "gabdt product with sum = " << (_gabdt_list_[index] * _sum_gabdt_list_[index - 1]) << std::endl
					<< "gabdt product with self = " << (_gabdt_list_[index] * _gabdt_list_[index]) << std::endl;


			std::cerr.flush();
		}
		throw std::runtime_error("ERROR: Could not calculate delta_rho in stripping_orbit_segment::_delta_rho.\n");
	}
	return max( result, 0. );
}

const double brgastro::stripping_orbit_segment::_step_length_factor( CONST_BRG_VELOCITY_REF  v, CONST_BRG_DISTANCE_REF  r ) const
{
	/* Method to determine the appropriate factor for the step length.
	 * It uses the host's virial velocity and radius so that parts of
	 * the orbit with either high velocity or low distance from the host
	 * will take small steps (and the inverse as well).
	 */
	BRG_VELOCITY v_0;
	BRG_DISTANCE r_0;

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

	if(isbad(step_length_factor_r))
	{
		std::cerr << "WARNING: Bad step_length_factor_r: " << step_length_factor_r
				<< "r = " << r
				<< "r_0 = " << r_0 << std::endl;
		if(isgood(step_length_factor_v)) return step_length_factor_v;
	}

	if(isbad(step_length_factor_v))
	{
		std::cerr << "WARNING: Bad step_length_factor_v: " << step_length_factor_r
				<< "v = " << r
				<< "v_0 = " << r_0 << std::endl;
		if(isgood(step_length_factor_r))
			return step_length_factor_v;
		else
			return 1;
	}

	// Use the smaller of the two calculated factors
	return min(step_length_factor_v,step_length_factor_r);
}

const BRG_DISTANCE brgastro::stripping_orbit_segment::_rvir(
		const int index ) const
{
	return _current_satellite_ptr_->rvir();
}
const int brgastro::stripping_orbit_segment::_pass_interpolation_type() const
{
	brgastro::interpolator::allowed_interpolation_type type_to_pass = brgastro::interpolator::SPLINE;
	if(_interpolation_type_ == brgastro::stripping_orbit::UNSET)
	{
		if(length() * std::sqrt(step_factor_min()) > _spline_resolution_) type_to_pass = brgastro::interpolator::LINEAR;
		if(length() * step_factor_min() > _spline_resolution_) type_to_pass = brgastro::interpolator::LOWER;
	}
	else
	{
		if(_interpolation_type_ == brgastro::stripping_orbit::LINEAR)
			type_to_pass = brgastro::interpolator::LINEAR;
		if(_interpolation_type_ == brgastro::stripping_orbit::SPLINE)
			type_to_pass = brgastro::interpolator::SPLINE;
		if(_interpolation_type_ == brgastro::stripping_orbit::UPPER)
			type_to_pass = brgastro::interpolator::UPPER;
		if(_interpolation_type_ == brgastro::stripping_orbit::LOWER)
			type_to_pass = brgastro::interpolator::LOWER;
	}

	int err_code = 0;

	_x_spline_.set_interpolation_type(type_to_pass);
	_y_spline_.set_interpolation_type(type_to_pass);
	_z_spline_.set_interpolation_type(type_to_pass);
	err_code += _vx_spline_.set_interpolation_type(type_to_pass);
	err_code += _vy_spline_.set_interpolation_type(type_to_pass);
	err_code += _vz_spline_.set_interpolation_type(type_to_pass);
	for(unsigned int i=0; i<_host_parameter_splines_.size(); i++)
	{
		_host_parameter_splines_[i].set_interpolation_type(type_to_pass);
	}
	if(err_code != 0)
		return LOWER_LEVEL_ERROR + err_code;
	return 0;
}

const int brgastro::stripping_orbit_segment::get_final_m_ret(
BRG_MASS & m_ret, const bool silent ) const
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
	m_ret = _init_satellite_ptr_->mtot()*_m_ret_list_.back();
	return 0;
}
const int brgastro::stripping_orbit_segment::get_final_frac_m_ret( double & frac_m_ret,
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
	frac_m_ret = _m_ret_list_.back() / safe_d(_m_ret_list_.front());
	return 0;
}
const int brgastro::stripping_orbit_segment::get_m_ret_points(
		std::vector< std::pair<double,double> > & m_ret_points,
		const bool silent ) const
{
	if ( (!_calculated_) || (!_record_full_data_) )
	{
		set_record_full_data(true);
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}

	m_ret_points.clear();
	for(size_t i = 0; i<_phase_output_list_.size(); i++)
	{
		m_ret_points.push_back( std::make_pair(
				_phase_output_list_[i].t,
				_m_ret_list_.at(i)));
	}

	return 0;
}

const int brgastro::stripping_orbit_segment::get_final_m_vir_ret(
BRG_MASS & m_ret, const bool silent ) const
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
	m_ret = _init_satellite_ptr_->mvir()*_m_vir_ret_list_.back();
	return 0;
}
const int brgastro::stripping_orbit_segment::get_final_frac_m_vir_ret( double & frac_m_ret,
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
	frac_m_ret = _m_vir_ret_list_.back() / safe_d(_m_vir_ret_list_.front());
	return 0;
}
const int brgastro::stripping_orbit_segment::get_m_vir_ret_points(
		std::vector< std::pair<double,double> > & m_ret_points,
		const bool silent ) const
{
	if ( (!_calculated_) || (!_record_full_data_) )
	{
		set_record_full_data(true);
		if ( calc() )
			return UNSPECIFIED_ERROR;
	}
	if( _bad_result_ )
	{
		return UNSPECIFIED_ERROR;
	}

	m_ret_points.clear();
	for(size_t i = 0; i<_phase_output_list_.size(); i++)
	{
		m_ret_points.push_back( std::make_pair(
				_phase_output_list_[i].t,
				_m_vir_ret_list_.at(i)));
	}

	return 0;
}

const int brgastro::stripping_orbit_segment::get_final_sum_deltarho(
long double & final_sum_deltarho, const bool silent ) const
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

const int brgastro::stripping_orbit_segment::get_final_sum_deltarho(
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
	final_sum_deltarho = (double)_sum_delta_rho_list_.back();
	return 0;
}

const int brgastro::stripping_orbit_segment::get_final_sum_gabdt(
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

const int brgastro::stripping_orbit_segment::clone_final_satellite(
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
				final_satellite_clone = new brgastro::tNFW_profile;
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

const int brgastro::stripping_orbit_segment::clone_final_host(
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
				final_host_clone = new brgastro::tNFW_profile;
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

const BRG_MASS brgastro::stripping_orbit_segment::final_m_ret() const
{
	BRG_MASS result = -1;

	if ( get_final_m_ret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_m_ret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_m_ret.\n");
	}
	return result;
}
const double brgastro::stripping_orbit_segment::final_frac_m_ret() const
{
	double result = -1;

	if ( get_final_frac_m_ret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_frac_m_ret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_frac_m_ret.\n");
	}
	return result;
}
const std::vector< std::pair<double,double> > brgastro::stripping_orbit_segment::m_ret_points() const
{
	std::vector< std::pair<double,double> > result;
	if( get_m_ret_points(result) )
	{
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::m_ret_points.\n");
	}

	return result;
}
const BRG_MASS brgastro::stripping_orbit_segment::final_m_vir_ret() const
{
	BRG_MASS result = -1;

	if ( get_final_m_vir_ret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_m_vir_ret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_m_vir_ret.\n");
	}
	return result;
}
const double brgastro::stripping_orbit_segment::final_frac_m_vir_ret() const
{
	double result = -1;

	if ( get_final_frac_m_vir_ret( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_frac_m_vir_ret.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_frac_m_vir_ret.\n");
	}
	return result;
}
const std::vector< std::pair<double,double> > brgastro::stripping_orbit_segment::m_vir_ret_points() const
{
	std::vector< std::pair<double,double> > result;
	if( get_m_vir_ret_points(result) )
	{
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::m_vir_ret_points.\n");
	}

	return result;
}

const long double brgastro::stripping_orbit_segment::final_sum_deltarho() const
{
	long double result;

	if ( get_final_sum_deltarho( result ) )
	{
		std::cerr << "ERROR: Could not calculate in stripping_orbit_segment::final_sum_deltarho.\n";
		std::cerr.flush();
		throw std::runtime_error("ERROR: Could not calculate in stripping_orbit_segment::final_sum_deltarho.\n");
		result = -1;
	}
	return result;
}

const brgastro::gabdt brgastro::stripping_orbit_segment::final_sum_gabdt() const
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

const bool & brgastro::stripping_orbit_segment::likely_disrupted() const
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

const brgastro::density_profile * brgastro::stripping_orbit_segment::final_satellite() const
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
			return NULL;
		}
	}
	if( _bad_result_ || _likely_disrupted_ )
		return _init_satellite_ptr_; // Sanest option in this case
	else
		return _current_satellite_ptr_;
}

const brgastro::density_profile * brgastro::stripping_orbit_segment::final_host() const
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
			return NULL;
		}
	}
	if( _bad_result_ || _likely_disrupted_ )
		return _init_host_ptr_; // Sanest option in this case
	else
		return _current_host_ptr_;
}

#endif

#endif // end brgastro::stripping_orbit_segment class function definitions
