/**********************************************************************\
  @file SALTSA_orbit.h

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

// body file: SALTSA_orbit.cpp

#ifndef __SALTSA_H_INCLUDED__
#include "SALTSA.h"
#else

#ifndef __SALTSA_ORBIT_H_INCLUDED__
#define __SALTSA_ORBIT_H_INCLUDED__

#include <vector>
#include <iostream>

#include "SALTSA_functor.hpp"
#include "SALTSA_phase.hpp"
#include "SALTSA.h"

namespace SALTSA {

/** Class Forward Declarations **/
#if (1)
// These declarations allow the various classes to point to each other without worry about
// which order they're declared in. (Order does still matter for classes containing other
// classes, though.)
//
// Classes are documented further in the definitions
class stripping_orbit;
class stripping_orbit_segment;

class interpolator_derivative;
class gabdt;

class interpolator_functor;
class interpolator_derivative_functor;
class interpolator_derivative_weight_functor;
class solve_rt_it_functor;
class solve_rt_grid_functor;
class gabdt_functor;

#endif // end Class forward declarations

/** Class Definitions **/
#if (1)

class interpolator_functor: public functor< double >
{
	/************************************************************
	 interpolator_functor
	 ---------------

	 Child of functor

	 This class is used to provide a functor * for getting
	 points along a spline (which is used here to differentiate the
	 spline).

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	const SALTSA::interpolator *_interpolator_ptr_;
	bool _interpolator_ptr_set_up_;

public:

	// Swap functions
	void swap(interpolator_functor & other);
	friend void swap(interpolator_functor & same, interpolator_functor & other) {same.swap(other);}

	// Constructors
	interpolator_functor();
	interpolator_functor(const interpolator_functor& other);
	interpolator_functor(const SALTSA::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_functor()
	{
	}

	// Operator=
	interpolator_functor& operator=(interpolator_functor other);

	// Set functions
	const int set_interpolator_ptr( const SALTSA::interpolator *new_interpolator_ptr );

	// Function method
	const int operator()( const double & in_param,
			double & out_param, const bool silent = false ) const;

};
// class interpolator_functor

class interpolator_derivative_functor: public functor< double >
{
	/************************************************************
	 interpolator_derivative_functor
	 --------------------------

	 Child of functor

	 This class is used to provide a functor * for getting
	 the derivative of an interpolator at a given point.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	interpolator_functor _interpolator_functor_;
	bool _interpolator_functor_set_up_;

public:

	// Swap functions
	void swap(interpolator_derivative_functor& other);
	friend void swap(interpolator_derivative_functor& same, interpolator_derivative_functor& other)
		{same.swap(other);}

	// Constructors
	interpolator_derivative_functor();
	interpolator_derivative_functor(const interpolator_derivative_functor& other);
	interpolator_derivative_functor( SALTSA::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_derivative_functor()
	{
	}

	// Operator=
	interpolator_derivative_functor& operator=(interpolator_derivative_functor other);

	// Set functions
	const int set_interpolator_ptr( const SALTSA::interpolator *new_interpolator_ptr );

	// Function method
	const int operator()( const double & in_param,
			double & out_param, const bool silent = false ) const;

};
// class interpolator_derivative_functor

class interpolator_derivative_weight_functor: public functor< double >
{
	/************************************************************
	 interpolator_derivative_weight_functor
	 ---------------------------------

	 Child of functor

	 This class is used to provide a functor * for getting
	 the weight of various points in the smoothing kernel for
	 calculating the derivative of an interpolator.

	 Use of this class is handled by the interpolator_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	double _sample_scale_, _sample_max_width_;
	double _t_min_, _t_max_, _centre_point_;

public:

	// Swap functions
	void swap(interpolator_derivative_weight_functor &other);
	friend void swap(interpolator_derivative_weight_functor &same,
			interpolator_derivative_weight_functor &other)
				{same.swap(other);}

	// Constructors
	interpolator_derivative_weight_functor();
	interpolator_derivative_weight_functor(const interpolator_derivative_weight_functor &other);

	// Destructor
	virtual ~interpolator_derivative_weight_functor()
	{
	}

	// Operator=
	interpolator_derivative_weight_functor & operator=(interpolator_derivative_weight_functor other);

	// Set functions
	const int set_sample_scale( double new_sample_scale );

	const int set_sample_max_width( double new_sample_max_width );

	const int set_center_point( double new_center_point );

	const int set_t_min( double new_t_min );

	const int set_t_max( double new_t_max );

	// Function method
	const int operator()( const double & in_param,
			double & out_param, const bool silent = false ) const;
};
// class interpolator_derivative_sample_functor

class interpolator_derivative
{
	/************************************************************
	 interpolator_derivative
	 -----------------

	 This class operates like an interpolator with added
	 features. The same functions can be used as with the basic
	 interpolator class, but this class also has the ability to point
	 to a different interpolator, which this is intended to be
	 the derivative of. Then, "unknown" domain points can be
	 passed to this class, and it will calculate the derivative of
	 the other interpolator at those points to help fill in gaps.

	 The points passed to this can all be "known" (domain and range
	 passed), all "unknown" (only domain passed, range calculated),
	 or a mix of the two.

	 Unknown points are calculated using a smoothing kernel to help
	 handle noise in the pointed-to interpolator. This must be adjusted
	 by hand based on how noisy the interpolator's points are for optimal
	 results. The noisier it is, the wider the kernel should be.
	 Use the set_sample_scale(...) and set_sample_max_width(...)
	 functions to adjust the kernel size (the scale is the sigma of
	 a Gaussian, the max_width is the limits of integration). Both
	 take in values as representing a fraction of t_max - t_min.

	 \************************************************************/
private:
	SALTSA::interpolator *_interpolator_ptr_;
	mutable SALTSA::interpolator _known_interpolator_, _estimated_interpolator_;
	bool _interpolator_ptr_set_up_;
	mutable bool _calculated_;

	interpolator_functor _interpolator_func_;

	std::vector< double > _unknown_t_list_;

	static double _default_sample_scale_, _default_sample_max_width_,
			_default_sample_precision_;
	double _sample_scale_, _sample_max_width_, _sample_precision_;
	mutable double _t_min_, _t_max_;

	SALTSA::interpolator::allowed_interpolation_type _interpolation_type_;

public:

	// Swap functions
	void swap(interpolator_derivative &other);
	friend void swap(interpolator_derivative &same, interpolator_derivative &other) {same.swap(other);}

	// Constructors
	interpolator_derivative();
	interpolator_derivative( const interpolator_derivative &other);
	interpolator_derivative( SALTSA::interpolator *init_interpolator_ptr );

	// Destructors
	virtual ~interpolator_derivative()
	{
	}

	// Operator=
	interpolator_derivative & operator=(interpolator_derivative other);

	// Set functions
	const int set_spline_ptr( SALTSA::interpolator *new_interpolator_ptr );
	const int clear_spline_ptr();
	const int set_default_sample_scale( double new_default_sample_scale );
	const int set_default_sample_max_width(
			double new_default_sample_max_width );
	const int set_sample_scale( double new_sample_scale );
	const int set_sample_max_width( double new_sample_max_width );
	const int reset_sample_scale(); // Sets it to default
	const int reset_sample_max_width(); // Sets it to default
	const int set_interpolation_type(
			SALTSA::interpolator::allowed_interpolation_type new_interpolation_type);
	const int reset_interpolation_type();

	// Functions for adding/clearing points
	const int add_point( const double t, const double x );
	const int add_unknown_point( const double t );

	const int clear_known_points();
	const int clear_unknown_points();
	const int clear_points(); // Clears all points

	// Full clear function
	const int clear();

	// Get functions
	const double operator()( double xval ) const;
};
// class interpolator_derivative

class gabdt
{
	/************************************************************
	 gabdt
	 -----

	 This class represents a useful physical construct, of how
	 much particles in halos have been disrupted by tidal shocking.

	 See the g_ab object from equation 10 of Taylor and Babul
	 (2001). This represents that multiplied by the timestep.

	 \************************************************************/

private:
	mutable bool _is_cached_;
	const density_profile *_host_ptr_;

	double _x_, _y_, _z_, _r_;
	double _dt_;
	mutable std::vector< std::vector< long double > > _dv_;

public:

	// Swap functions
	void swap(gabdt &other);
	friend void swap(gabdt &same, gabdt &other) {same.swap(other);}

	// Constructors
	gabdt();
	gabdt( const gabdt & other_gabdt );
	gabdt( const density_profile *init_host, const double init_x,
			const double init_y, const double init_z,
			const double init_dt );

	// Destructor
	virtual ~gabdt();

	// Full clear function
	const int clear();

	// Set functions
	const int set( const SALTSA::density_profile *new_host_ptr,
			const double new_x, const double new_y,
			const double new_z, const double new_dt );
	const int set_pos( const double new_x, const double new_y,
			const double new_z );
	const int set_dt( const double dt );
	const int set_host_ptr( const density_profile *new_host_ptr );
	const int override_zero();

	// Calculation function
	const int calc_dv( const bool silent = false ) const;

	// Get functions
	const density_profile * host() const;
	const double x() const;
	const double y() const;
	const double z() const;
	const double r() const;
	const std::vector< std::vector< long double > > dv() const; // involves calculation if necessary
	const long double dv( const int x_i, const int y_i ) const; // involves calculation if necessary

	// Operator overloading
	const double operator*( const gabdt & other_gabdt ) const; // Dot-product(ish) operator

	gabdt & operator=( gabdt other_gabdt ); // Assignment

	gabdt & operator+=( const gabdt & other_gabdt ); // Addition
	gabdt operator+( const gabdt & other_gabdt ) const;

	gabdt & operator*=( const double scale_fraction ); // Multiplication by a double
	gabdt operator*( const double scale_fraction ) const;

};

class gabdt_functor: public functor< std::vector< double > >
{
	/************************************************************
	 gabdt_functor
	 --------------

	 Child of functor

	 This class provides a functor * for getting the 3-D
	 acceleration within a halo.

	 \************************************************************/
public:

	// Swap functions
	void swap(gabdt_functor &other);
	friend void swap(gabdt_functor &same, gabdt_functor &other) {same.swap(other);}

	// Constructors
	gabdt_functor();
	gabdt_functor(const gabdt_functor &other);

	// Destructor
	virtual ~gabdt_functor()
	{
	}

	// Operator=
	gabdt_functor & operator=(gabdt_functor other);

	// Host accessor
	const density_profile *host_ptr;

	// Function method
	const int operator()( const std::vector< double > & in_params,
			std::vector< double > & out_params,
			const bool silent = false ) const;

}; // class gabdt_function


class stripping_orbit_segment
{
	/************************************************************
	 stripping_orbit_segment
	 -----------------------

	 An individual segment of a stripping_orbit which has been
	 split up by discontinuity times. This class can actually be
	 used just like a stripping_orbit with the exception of
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
	double _t_min_natural_value_, _t_max_natural_value_, _t_min_override_val_, _t_max_override_val_;
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
	mutable SALTSA::interpolator _x_spline_, _y_spline_, _z_spline_,
			_test_mass_spline_;
	mutable interpolator_derivative _vx_spline_, _vy_spline_, _vz_spline_;
	mutable std::vector< SALTSA::interpolator > _host_parameter_splines_;

	// Vectors for output data
	mutable std::vector< double > _delta_rho_list_,
			_x_data_, _y_data_, _z_data_, _vx_data_, _vy_data_, _vz_data_,
			_rt_list_, _rt_ratio_list_;
	mutable std::vector< long double > _sum_delta_rho_list_, _m_ret_list_;
	mutable std::vector< double > _m_vir_ret_list_;
	mutable std::vector< std::vector< double > > _satellite_parameter_data_; // Keeps track of satellite's parameters (ie. mass, tau)
	mutable std::vector< std::vector< double > > _host_parameter_data_;

	// Other maintained vectors for calculations
	mutable std::vector< gabdt > _gabdt_list_;
	mutable std::vector< gabdt > _sum_gabdt_list_;
	mutable std::vector< SALTSA::phase > _phase_list_, _phase_output_list_;

	// Output-modifying parameters
	int _num_parameters_;
	mutable std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;

	// Private functions
	const int _init();
	const int _reserve( int n, const bool silent = false ) const;
	const double _delta_rho( const int index, const double x,
			const double t_step, const bool silent = false ) const;
	const double _step_length_factor( const double  v, const double  r ) const;
	const double _rvir( const int index = 0 ) const;
	const int _pass_interpolation_type() const;

	// Calculation assistance functions
	const double _tidal_strip_retained( const density_profile *host,
			const density_profile *satellite, const double r,
			const double vr, const double vt,
			const double time_step, const long double &sum_delta_rho = 0 ) const;
	const double _get_rt( const density_profile *host,
			const density_profile *satellite, const double r,
			const double vr, const double vt,
			const double time_step, const long double &sum_delta_rho,
			const bool silent = false ) const;
#endif

public:

	// Swap functions
	void swap(stripping_orbit_segment &other);
	friend void swap(stripping_orbit_segment &same, stripping_orbit_segment &other) {same.swap(other);}

	// Constructors, destructors, and related operations
	stripping_orbit_segment();
	stripping_orbit_segment(
			const stripping_orbit_segment &other );
	stripping_orbit_segment & operator=(
			stripping_orbit_segment other );
	stripping_orbit_segment( const density_profile *host,
			const density_profile *satellite,
			const int init_resolution = 200 );
	virtual ~stripping_orbit_segment();
	stripping_orbit_segment *stripping_orbit_spline_clone() const;

	// Functions to add points to the splines
	const int add_point( const double x, const double y,
			const double z, const double t,
			double new_test_mass = 1 );
	const int add_point( const double x, const double y,
			const double z, const double vx,
			const double vy, const double vz, const double t,
			double new_test_mass = 1 );
	const int add_x_point( const double x, const double t );
	const int add_y_point( const double y, const double t );
	const int add_z_point( const double z, const double t );
	const int add_vx_point( const double vx, const double t );
	const int add_vy_point( const double vy, const double t );
	const int add_vz_point( const double vz, const double t );
	const int add_unknown_vx_point( const double t );
	const int add_unknown_vy_point( const double t );
	const int add_unknown_vz_point( const double t );
	const int add_test_mass_point( const double test_mass, const double t );
	const int add_host_parameter_point( const unsigned int num_parameters,
			const std::vector< double > &parameters, const double t,
			const bool silent = false );


	// Setting integration parameters
#if(1)
	const int set_resolution( const int new_spline_resolution,
			const bool silent=false );
	const int set_interpolation_type( const stripping_orbit::allowed_interpolation_type new_type,
			const bool silent=false );
	const int set_v_0( const double new_v_0,
			const bool silent=false );
	const int set_r_0( const double new_r_0,
			const bool silent=false );
	const int set_step_length_power( const double new_step_length_power,
			const bool silent=false );
	const int set_step_factor_max( const double new_step_factor_max,
			const bool silent=false );
	const int set_step_factor_min( const double new_step_factor_min,
			const bool silent=false );
#endif

	// Setting tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	const int set_tidal_stripping_amplification( const double new_tidal_stripping_amplification,
			const bool silent=false );
	const int set_tidal_stripping_deceleration( const double new_tidal_stripping_deceleration,
			const bool silent=false );
	const int set_tidal_stripping_radialness( const double new_tidal_stripping_radialness,
			const bool silent=false );
	const int set_tidal_shocking_amplification( const double new_tidal_shocking_amplification,
			const bool silent=false );
	const int set_tidal_shocking_persistance( const double new_tidal_shocking_persistance,
			const bool silent=false );
	const int set_tidal_shocking_power( const double new_tidal_shocking_power,
			const bool silent=false );
#endif


	// Set initial/global parameters
	const int set_tNFW_init_satellite( const double new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 );
	const int set_tNFW_host( const double new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 );
	const int set_t_min( const double new_t_min );
	const int set_t_max( const double new_t_max );
	const int reset_t_min();
	const int reset_t_max();
	const int set_init_sum_deltarho( const long double &new_init_sum_deltarho );
	const int set_init_sum_gabdt( const gabdt &new_init_gabdt );
	const int set_init_satellite( const density_profile *new_init_satellite );
	const int set_init_host( const density_profile *new_init_host );

	// Clear orbit data or initial parameters
	const int clear_points();
	const int clear_host_parameter_points();
	const int clear_init_sum_deltarho();
	const int clear_init_sum_gabdt();
	const int clear_init_satellite();
	const int clear_init_host();

	// Global clearing functions
	const int clear();
	const int clear_calcs() const; // Only clears calculations

	// Functions for determining how calc() will be called
	const int set_record_full_data( const bool new_record_full_data ) const;

	// Function to calculate stripping
	const int calc( const bool silent = false ) const;

	// Output-modifying functions
	const int set_satellite_output_parameters(
			const unsigned int num_parameters,
			const std::vector< bool > &satellite_output_parameters );
	const int set_satellite_parameter_unitconvs(
			const unsigned int num_parameters,
			const std::vector< double > &satellite_unitconvs );
	const int set_host_output_parameters( const unsigned int num_parameters,
			const std::vector< bool > &host_output_parameters );
	const int set_host_parameter_unitconvs( const unsigned int num_parameters,
			const std::vector< double > &host_unitconvs );
	const int clear_satellite_output_parameters();
	const int clear_satellite_parameter_unitconvs();
	const int clear_host_output_parameters();
	const int clear_host_parameter_unitconvs();

	// Print output function
	const int print_full_data( std::ostream *out, const bool include_header =
			true, const double m_ret_multiplier = 1, const double m_vir_ret_multiplier = 1,
			const bool silent = false ) const;

	// Get current state of object
	const unsigned int length() const;

	// Accessors
#if(1)

	// Integration parameters
#if(1)
	const int & spline_resolution() const {return _spline_resolution_;}
	const stripping_orbit::allowed_interpolation_type & interpolation_type()
		{return _interpolation_type_;}
	const double  v_0() const {return _v_0_;}
	const double  r_0() const {return _r_0_;}
	const double  step_length_power() const {return _step_length_power_;}
	const double  step_factor_max() const {return _step_factor_max_;}
	const double  step_factor_min() const {return _step_factor_min_;}
#endif

	// Tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	const double  tidal_stripping_amplification() const {return _tidal_stripping_amplification_;}
	const double  tidal_stripping_deceleration() const {return _tidal_stripping_deceleration_;}
	const double  tidal_stripping_radialness() const {return _tidal_stripping_radialness_;}
	const double  tidal_shocking_amplification() const {return _tidal_shocking_amplification_;}
	const double  tidal_shocking_persistance() const {return _tidal_shocking_persistance_;}
	const double  tidal_shocking_power() const {return _tidal_shocking_power_;}
#endif

	const bool & calculated() const {return _calculated_;};
	const bool & bad_result() const {return _bad_result_;};
	const density_profile * init_satellite_ptr() const {return _init_satellite_ptr_;};
	const density_profile * init_host_ptr() const {return _init_host_ptr_;};
	const double  t_min_natural_value() const {return _t_min_natural_value_;};
#endif

	// Get final data (returns 1 on error)
	const int get_final_m_ret( double & m_ret,
			const bool silent = false ) const;
	const int get_final_frac_m_ret( double & final_frac_m_ret,
			const bool silent = false ) const;
	const int get_final_m_vir_ret( double & m_vir_ret,
			const bool silent = false ) const;
	const int get_final_frac_m_vir_ret( double & final_frac_m_vir_ret,
			const bool silent = false ) const;
	const int get_final_sum_deltarho( long double & final_sum_deltarho,
			const bool silent = false ) const;
	const int get_final_sum_deltarho( double & final_sum_deltarho,
			const bool silent = false ) const;
	const int get_m_ret_points( std::vector< std::pair<double,double> > & m_ret_points,
			const bool silent = false ) const;
	const int get_m_vir_ret_points( std::vector< std::pair<double,double> > & m_vir_ret_points,
			const bool silent = false ) const;
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt, const bool silent =
			false ) const;
	const int clone_final_satellite( density_profile * & final_satellite_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!
	const int clone_final_host( density_profile * & final_host_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!

	// Get final data (throws exception on error)
	const double final_m_ret() const;
	const double final_frac_m_ret() const;
	const double final_m_vir_ret() const;
	const double final_frac_m_vir_ret() const;
	const long double final_sum_deltarho() const;
	const std::vector< std::pair<double,double> > m_ret_points() const;
	const std::vector< std::pair<double,double> > m_vir_ret_points() const;
	const gabdt final_sum_gabdt() const;
	const bool & likely_disrupted() const;
	const density_profile * final_satellite() const; // Creates a clone. Make sure to delete!
	const density_profile * final_host() const; // Creates a clone. Make sure to delete!

};
// class stripping_orbit_segment

class solve_rt_it_functor: public functor< double > // Always uses one param, returns new value for that parameter for iteration.
{
	/************************************************************
	 solve_rt_it_functor
	 --------------------

	 Child of functor

	 Provides a functor * to be used to iteratively solve
	 for tidal radius.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	long double sum_delta_rho, Daccel, omega;
	const int operator()( const double & in_param,
			double & out_param, const bool silent = false ) const;
	solve_rt_it_functor( const double init_omega,
			const density_profile *init_satellite, const double init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_it_functor();
};

class solve_rt_grid_functor: public functor< double >
{
	/************************************************************
	 solve_rt_grid_functor
	 ----------------------

	 Child of functor

	 Provides a functor * to be used with the grid solver.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	long double sum_delta_rho, Daccel, omega;
	const int operator()( const double & in_param,
			double & out_param, const bool silent = false ) const;
	solve_rt_grid_functor( const double init_omega,
			const density_profile *init_satellite, const double init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_grid_functor();
};

#endif // end class definitions

} // end namespace SALTSA

#endif // __SALTSA_ORBIT_H_INCLUDED__

#endif // __SALTSA_H_INCLUDED__
