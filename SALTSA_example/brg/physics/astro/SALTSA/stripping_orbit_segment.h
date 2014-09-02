/**********************************************************************\
  @file stripping_orbit_segment.h

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

// body file: brg/physics/astro/SALTSA/stripping_orbit_segment.cpp

#ifndef _STRIPPING_ORBIT_H_INCLUDED_
#include "brg/physics/astro/SALTSA/stripping_orbit.h"
#else

#ifndef _STRIPPING_ORBIT_SEGMENT_H_INCLUDED_
#define _STRIPPING_ORBIT_SEGMENT_H_INCLUDED_

#include <vector>

#include "brg/global.h"

#include "brg/math/interpolator/interpolator.h"
#include "brg/math/interpolator/interpolator_derivative.h"

#include "brg/physics/astro/SALTSA/gabdt.h"
#include "brg/physics/units/unit_obj.h"

namespace brgastro {

class stripping_orbit_segment
{
	/************************************************************
	 stripping_orbit_segment
	 -----------------------

	 An individual segment of a stripping_orbit which has been
	 split up by discontinuity times. This class can actually be
	 used much like a stripping_orbit with the exception of
	 adding discontinuity times to it (this will be marginally
	 faster, but not really noticeable).

	 \************************************************************/
private:
#if(1)

	// Integration parameters
#if(1)
	int _spline_resolution_;

	// Interpolation method
	stripping_orbit::allowed_interpolation_type _interpolation_type_;

	// Variable step length tweaking: Time step length is proportional to (v_0/v)^(step_length_power)
	// This gives smaller steps when the satellite is moving faster.
	// If you want to turn off adaptive step size, set step_length_power to 0
	// Alternatively, set step_length_power to 1 for even steps in position
	double _v_0_; // 400 km/s
	double _r_0_; // 400 kpc
	double _step_length_power_; // How strongly variable step length is implemented
	double _step_factor_max_; // Maximum allowed value of (v_0/v)^(step_length_power)
	double _step_factor_min_; // Minimum allowed value of (v_0/v)^(step_length_power)
#endif

	// Tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	double _tidal_stripping_amplification_; // Amplifies tidal stripping by this factor
	double _tidal_stripping_deceleration_; // If positive, increase tidal stripping near pericentre,
	  	  	  	  	  	  	  	  	  	   // if negative, decrease near pericentre
	double _tidal_stripping_radialness_; // How much tidal stripping depends on full velocity
	                                     // versus tangential velocity. Larger value of this
	                                     // increases stripping of more radial orbits preferentially
	double _tidal_shocking_amplification_; // Amplifies tidal heating by this factor
	double _tidal_shocking_persistance_; // How long shocking is active for
	double _tidal_shocking_power_; // Affects interplay of stripping and satellite halo profile
#endif

	// Initial parameters
	const density_profile *_init_host_ptr_, *_init_satellite_ptr_;
	mutable density_profile *_current_host_ptr_, *_current_satellite_ptr_;
	long double _init_sum_delta_rho_;
	gabdt _init_sum_gabdt_;

	// Global parameters
	BRG_TIME _t_min_natural_value_, _t_max_natural_value_, _t_min_override_val_, _t_max_override_val_;
	bool _override_t_min_, _override_t_max_;
	mutable bool _record_full_data_;
	bool _host_loaded_, _satellite_loaded_;
	mutable bool _calculated_, _bad_result_, _current_satellite_in_use_,
			_current_host_in_use_, _likely_disrupted_;
	bool _evolving_host_;
	bool _using_private_init_host_, _using_private_init_satellite_;
	tNFW_profile _private_tNFW_init_host_, _private_tNFW_init_satellite_;

	// Splines for orbit data
	// Must be mutable since the spline class used doesn't allow conceptual constness for calculating
	mutable brgastro::interpolator _x_spline_, _y_spline_, _z_spline_,
			_test_mass_spline_;
	mutable interpolator_derivative _vx_spline_, _vy_spline_, _vz_spline_;
	mutable std::vector< brgastro::interpolator > _host_parameter_splines_;

	// Vectors for output data
	mutable std::vector< BRG_UNITS > _delta_rho_list_,
			_x_data_, _y_data_, _z_data_, _vx_data_, _vy_data_, _vz_data_,
			_rt_list_, _rt_ratio_list_;
	mutable std::vector< long double > _sum_delta_rho_list_, _m_ret_list_;
	mutable std::vector< double >_m_vir_ret_list_;
	mutable std::vector< std::vector< BRG_UNITS > > _satellite_parameter_data_; // Keeps track of satellite's parameters (ie. mass, tau)
	mutable std::vector< std::vector< BRG_UNITS > > _host_parameter_data_;

	// Other maintained vectors for calculations
	mutable std::vector< gabdt > _gabdt_list_;
	mutable std::vector< gabdt > _sum_gabdt_list_;
	mutable std::vector< brgastro::phase > _phase_list_, _phase_output_list_;

	// Output-modifying parameters
	unsigned int _num_parameters_;
	mutable std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;

	// Private functions
	void _init();
	void _reserve( const unsigned int n, const bool silent = false ) const;
	BRG_UNITS _delta_rho( const int index, const double x,
			CONST_BRG_TIME_REF t_step, const bool silent = false ) const;
	const double _step_length_factor( CONST_BRG_VELOCITY_REF  v, CONST_BRG_DISTANCE_REF  r ) const;
	BRG_DISTANCE _rvir( const int index = 0 ) const;
	void _pass_interpolation_type() const;

	// Calculation assistance functions
	double _tidal_strip_retained( const density_profile *host,
			const density_profile *satellite, CONST_BRG_DISTANCE_REF r,
			CONST_BRG_VELOCITY_REF vr, CONST_BRG_VELOCITY_REF vt,
			CONST_BRG_TIME_REF time_step, const long double &sum_delta_rho = 0 ) const;
	BRG_DISTANCE _get_rt( const density_profile *host,
			const density_profile *satellite, CONST_BRG_DISTANCE_REF r,
			CONST_BRG_VELOCITY_REF vr, CONST_BRG_VELOCITY_REF vt,
			CONST_BRG_TIME_REF time_step, const long double &sum_delta_rho,
			const bool silent = false ) const;
#endif

public:

	// Swap function
	void swap(stripping_orbit_segment &other);

	// Constructors, destructors, and related operations
	stripping_orbit_segment();
	stripping_orbit_segment(
			const stripping_orbit_segment &other );
	stripping_orbit_segment & operator=(
			stripping_orbit_segment other );
	stripping_orbit_segment( const density_profile *host,
			const density_profile *satellite,
			const unsigned int init_resolution = 200 );
	virtual ~stripping_orbit_segment();
	stripping_orbit_segment *stripping_orbit_spline_clone() const;

	// Move constructor and move assignment (C++11 only)
#ifdef _BRG_USE_CPP_11_STD_
	stripping_orbit_segment(stripping_orbit_segment &&other)
	: stripping_orbit_segment()
	{
		swap(other);
	}
	stripping_orbit_segment & operator=(stripping_orbit_segment &&other)
	{
		swap(other);
		return *this;
	}
#endif

	// Functions to add points to the splines
	void add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
			double new_test_mass = 1 );
	void add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
			CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
			double new_test_mass = 1 );
	void add_x_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_TIME_REF t );
	void add_y_point( CONST_BRG_DISTANCE_REF y, CONST_BRG_TIME_REF t );
	void add_z_point( CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t );
	void add_vx_point( CONST_BRG_VELOCITY_REF vx, CONST_BRG_TIME_REF t );
	void add_vy_point( CONST_BRG_VELOCITY_REF vy, CONST_BRG_TIME_REF t );
	void add_vz_point( CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t );
	void add_unknown_vx_point( CONST_BRG_TIME_REF t );
	void add_unknown_vy_point( CONST_BRG_TIME_REF t );
	void add_unknown_vz_point( CONST_BRG_TIME_REF t );
	void add_test_mass_point( const double test_mass, CONST_BRG_TIME_REF t );
	void add_host_parameter_point(
			const std::vector< BRG_UNITS > &parameters, CONST_BRG_TIME_REF t,
			const bool silent = false );


	// Setting integration parameters
#if(1)
	void set_resolution( const unsigned int new_spline_resolution,
			const bool silent=false );
	void set_interpolation_type( const stripping_orbit::allowed_interpolation_type new_type,
			const bool silent=false );
	void set_v_0( const double new_v_0,
			const bool silent=false );
	void set_r_0( const double new_r_0,
			const bool silent=false );
	void set_step_length_power( const double new_step_length_power,
			const bool silent=false );
	void set_step_factor_max( const double new_step_factor_max,
			const bool silent=false );
	void set_step_factor_min( const double new_step_factor_min,
			const bool silent=false );
#endif

	// Setting tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	void set_tidal_stripping_amplification( const double new_tidal_stripping_amplification,
			const bool silent=false );
	void set_tidal_stripping_deceleration( const double new_tidal_stripping_deceleration,
			const bool silent=false );
	void set_tidal_stripping_radialness( const double new_tidal_stripping_radialness,
			const bool silent=false );
	void set_tidal_shocking_amplification( const double new_tidal_shocking_amplification,
			const bool silent=false );
	void set_tidal_shocking_persistance( const double new_tidal_shocking_persistance,
			const bool silent=false );
	void set_tidal_shocking_power( const double new_tidal_shocking_power,
			const bool silent=false );
#endif


	// Set initial/global parameters
	void set_tNFW_init_satellite( CONST_BRG_MASS_REF new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 );
	void set_tNFW_host( CONST_BRG_MASS_REF new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 );
	void set_t_min( CONST_BRG_TIME_REF new_t_min );
	void set_t_max( CONST_BRG_TIME_REF new_t_max );
	void reset_t_min();
	void reset_t_max();
	void set_init_sum_deltarho( const long double &new_init_sum_deltarho );
	void set_init_sum_gabdt( const gabdt &new_init_gabdt );
	void set_init_satellite( const density_profile *new_init_satellite );
	void set_init_host( const density_profile *new_init_host );

	// Clear orbit data or initial parameters
	void clear_points();
	void clear_host_parameter_points();
	void clear_init_sum_deltarho();
	void clear_init_sum_gabdt();
	void clear_init_satellite();
	void clear_init_host();

	// Global clearing functions
	void clear();
	void clear_calcs() const; // Only clears calculations

	// Functions for determining how calc() will be called
	void set_record_full_data( const bool new_record_full_data ) const;

	// Function to calculate stripping
	void calc( const bool silent = false ) const;

	// Output-modifying functions
	void set_satellite_output_parameters(
			const std::vector< bool > &satellite_output_parameters );
	void set_satellite_parameter_unitconvs(
			const std::vector< double > &satellite_unitconvs );
	void set_host_output_parameters(
			const std::vector< bool > &host_output_parameters );
	void set_host_parameter_unitconvs(
			const std::vector< double > &host_unitconvs );
	void clear_satellite_output_parameters();
	void clear_satellite_parameter_unitconvs();
	void clear_host_output_parameters();
	void clear_host_parameter_unitconvs();

	// Print output function
	void print_full_data( std::ostream *out, const bool include_header =
			true, const double m_ret_multiplier = 1, const double m_vir_ret_multiplier = 1,
			const bool silent = false ) const;

	// Get current state of object
	unsigned int length() const;

	// Accessors
#if(1)

	// Integration parameters
#if(1)
	unsigned int spline_resolution() const {return _spline_resolution_;}
	stripping_orbit::allowed_interpolation_type interpolation_type()
		{return _interpolation_type_;}
	BRG_VELOCITY v_0() const {return _v_0_;}
	BRG_DISTANCE r_0() const {return _r_0_;}
	double step_length_power() const {return _step_length_power_;}
	double step_factor_max() const {return _step_factor_max_;}
	double step_factor_min() const {return _step_factor_min_;}
#endif

	// Tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	double tidal_stripping_amplification() const {return _tidal_stripping_amplification_;}
	double tidal_stripping_deceleration() const {return _tidal_stripping_deceleration_;}
	double tidal_stripping_radialness() const {return _tidal_stripping_radialness_;}
	double tidal_shocking_amplification() const {return _tidal_shocking_amplification_;}
	double tidal_shocking_persistance() const {return _tidal_shocking_persistance_;}
	double tidal_shocking_power() const {return _tidal_shocking_power_;}
#endif

	bool calculated() const {return _calculated_;};
	bool bad_result() const {return _bad_result_;};
	const density_profile * init_satellite_ptr() const {return _init_satellite_ptr_;};
	const density_profile * init_host_ptr() const {return _init_host_ptr_;};
	BRG_TIME t_min_natural_value() const {return _t_min_natural_value_;};
#endif

	// Get final data (returns 1 on error)
	int get_final_m_ret( BRG_MASS & m_ret,
			const bool silent = false ) const;
	int get_final_frac_m_ret( double & final_frac_m_ret,
			const bool silent = false ) const;
	int get_final_m_vir_ret( BRG_MASS & m_vir_ret,
			const bool silent = false ) const;
	int get_final_frac_m_vir_ret( double & final_frac_m_vir_ret,
			const bool silent = false ) const;
	int get_final_sum_deltarho( long double & final_sum_deltarho,
			const bool silent = false ) const;
	int get_final_sum_deltarho( BRG_UNITS & final_sum_deltarho,
			const bool silent = false ) const;
	int get_m_ret_points( std::vector< std::pair<double,double> > & m_ret_points,
			const bool silent = false ) const;
	int get_m_vir_ret_points( std::vector< std::pair<double,double> > & m_vir_ret_points,
			const bool silent = false ) const;
	int get_final_sum_gabdt( gabdt & final_sum_gabdt, const bool silent =
			false ) const;
	int clone_final_satellite( density_profile * & final_satellite_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!
	int clone_final_host( density_profile * & final_host_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!

	// Get final data (throws exception on error)
	BRG_MASS final_m_ret() const;
	double final_frac_m_ret() const;
	BRG_MASS final_m_vir_ret() const;
	double final_frac_m_vir_ret() const;
	BRG_UNITS final_sum_deltarho() const;
	long double long_final_sum_deltarho() const;
	std::vector< std::pair<double,double> > m_ret_points() const;
	std::vector< std::pair<double,double> > m_vir_ret_points() const;
	gabdt final_sum_gabdt() const;
	bool likely_disrupted() const;
	const density_profile * final_satellite() const; // Creates a clone. Make sure to delete!
	const density_profile * final_host() const; // Creates a clone. Make sure to delete!

};
// class stripping_orbit_segment

} // end namespace brgastro



#endif /* _STRIPPING_ORBIT_SEGMENT_H_INCLUDED_ */

#endif // ifndef _STRIPPING_ORBIT_H_INCLUDED_
