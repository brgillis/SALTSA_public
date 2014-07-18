/**       @file SALTSA_orbit.h
 *
 *     Project: SALTSA_example
 *        Path: /SALTSA_example/SALTSA_orbit.h
 *
 *  Created on: 17 Jul 2014
 *      Author: brg
 */

#ifndef __SALTSA_ORBIT_H_INCLUDED__
#define __SALTSA_ORBIT_H_INCLUDED__

#include <vector>
#include <iostream>

#include "SALTSA_phase.hpp"

namespace SALTSA {

/** Class Forward Declarations **/
#if (1)
// These declarations allow the various classes to point to each other without worry about
// which order they're declared in. (Order does still matter for classes containing other
// classes, though.)
//
// Classes are documented further in the definitions
class stripping_orbit_segment;

class spline_derivative;
class gabdt;

class spline_function;
class spline_derivative_function;
class spline_derivative_weight_function;
class solve_rt_it_function;
class solve_rt_grid_function;
class gabdt_function;

#endif // end Class forward declarations

/** Class Definitions **/
#if (1)

class spline_function: public functor< double >
{
	/************************************************************
	 spline_function
	 ---------------

	 Child of function_class

	 This class is used to provide a function_class * for getting
	 points along a spline (which is used here to differentiate the
	 spline).

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	SALTSA::interpolator *_spline_ptr_;
	bool _spline_ptr_set_up_;

public:

	// Constructors
	spline_function();
	spline_function( SALTSA::interpolator *init_spline_ptr );

	// Destructor
	virtual ~spline_function()
	{
	}

	// Set functions
	const int set_spline_ptr( SALTSA::interpolator *new_spline_ptr );

	// Function method
	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;

};
// class spline_function

class spline_derivative_function: public functor< double >
{
	/************************************************************
	 spline_derivative_function
	 --------------------------

	 Child of function_class

	 This class is used to provide a function_class * for getting
	 the derivative of a spline at a given point.

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	spline_function _spline_function_;
	bool _spline_function_set_up_;

public:

	// Constructors
	spline_derivative_function();
	spline_derivative_function( SALTSA::interpolator *init_spline_ptr );

	// Destructor
	virtual ~spline_derivative_function()
	{
	}

	// Set functions
	const int set_spline_ptr( SALTSA::interpolator *new_spline_ptr );

	// Function method
	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;

};
// class spline_derivative_function

class spline_derivative_weight_function: public functor< double >
{
	/************************************************************
	 spline_derivative_weight_function
	 ---------------------------------

	 Child of function_class

	 This class is used to provide a function_class * for getting
	 the weight of various points in the smoothing kernel for
	 calculating the derivative of a spline.

	 Use of this class is handled by the spline_derivative class.
	 No need for the end-user to worry too much about how this
	 works.

	 \************************************************************/

private:

	double _sample_scale_, _sample_max_width_;
	double _t_min_, _t_max_, _centre_point_;

public:

	// Constructors
	spline_derivative_weight_function();

	// Destructor
	virtual ~spline_derivative_weight_function()
	{
	}

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
// class spline_derivative_sample_function

class spline_derivative
{
	/************************************************************
	 spline_derivative
	 -----------------

	 This class operates like a magnet::math::Spline with added
	 features. The same functions can be used as with the basic
	 Spline class, but this class also has the ability to point
	 to a different spline, which this spline is intended to be
	 the derivative of. Then, "unknown" domain points can be
	 passed to this class, and it will calculate the derivative of
	 the other spline at those points to help fill in gaps.

	 The points passed to this can all be "known" (domain and range
	 passed), all "unknown" (only domain passed, range calculated),
	 or a mix of the two.

	 Unknown points are calculated using a smoothing kernel to help
	 handle noise in the pointed-to spline. This must be adjusted
	 by hand based on how noisy the spline's points are for optimal
	 results. The noisier it is, the wider the kernel should be.
	 Use the set_sample_scale(...) and set_sample_max_width(...)
	 functions to adjust the kernel size (the scale is the sigma of
	 a Gaussian, the max_width is the limits of integration). Both
	 take in values as representing a fraction of t_max - t_min.

	 \************************************************************/
private:
	SALTSA::interpolator *_spline_ptr_;
	mutable SALTSA::interpolator _known_spline_, _estimated_spline_;
	bool _spline_ptr_set_up_;
	mutable bool _calculated_;

	spline_function _spline_func_;

	std::vector< double > _unknown_t_list_;

	static double _default_sample_scale_, _default_sample_max_width_,
			_default_sample_precision_;
	double _sample_scale_, _sample_max_width_, _sample_precision_;
	mutable double _t_min_, _t_max_;

public:
	// Constructors
	spline_derivative();
	spline_derivative( SALTSA::interpolator *init_spline_ptr );

	// Destructors
	virtual ~spline_derivative()
	{
	}

	// Set functions
	const int set_spline_ptr( SALTSA::interpolator *new_spline_ptr );
	const int clear_spline_ptr();
	const int set_default_sample_scale( double new_default_sample_scale );
	const int set_default_sample_max_width(
			double new_default_sample_max_width );
	const int set_sample_scale( double new_sample_scale );
	const int set_sample_max_width( double new_sample_max_width );
	const int reset_sample_scale(); // Sets it to default
	const int reset_sample_max_width(); // Sets it to default

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
// class spline_derivative

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

	double _x_, _y_, _z_, _r_;double _dt_;
	mutable std::vector< std::vector< double > > _dv_;

public:

	// Constructors
	gabdt();
	gabdt( const gabdt & other_gabdt );
	gabdt( const density_profile *init_host, const double &init_x,
			const double &init_y, const double &init_z,
			const double &init_dt );

	// Destructor
	virtual ~gabdt();

	// Full clear function
	const int clear();

	// Set functions
	const int set( const SALTSA::density_profile *new_host_ptr,
			const double &new_x, const double &new_y,
			const double &new_z, const double &new_dt );
	const int set_pos( const double &new_x, const double &new_y,
			const double &new_z );
	const int set_dt( const double &dt );
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
	const std::vector< std::vector< double > > dv() const; // involves calculation if necessary
	const double dv( const int x_i, const int y_i ) const; // involves calculation if necessary

	// Operator overloading
	const double operator*( const gabdt & other_gabdt ) const; // Dot-product(ish) operator

	gabdt & operator=( const gabdt & other_gabdt ); // Assignment

	gabdt & operator+=( const gabdt & other_gabdt ); // Addition
	gabdt operator+( const gabdt & other_gabdt ) const;

	gabdt & operator*=( const double scale_fraction ); // Multiplication by a double
	gabdt operator*( const double scale_fraction ) const;

};

class gabdt_function: public functor< std::vector< double > >
{
	/************************************************************
	 gabdt_function
	 --------------

	 Child of function_class

	 This class provides a function_class * for getting the 3-D
	 acceleration within a halo.

	 \************************************************************/
public:

	// Constructor
	gabdt_function();

	// Destructor
	virtual ~gabdt_function()
	{
	}

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
	double _tidal_shocking_amplification_; // Amplifies tidal heating by this factor
	double _tidal_shocking_persistance_; // How long shocking is active for
	double _tidal_shocking_power_; // Affects interplay of stripping and satellite halo profile
#endif

	// Initial parameters
	const density_profile *_init_host_ptr_, *_init_satellite_ptr_;
	mutable density_profile *_current_host_ptr_, *_current_satellite_ptr_;
	double _init_sum_delta_rho_;
	gabdt _init_sum_gabdt_;

	// Global parameters
	double _t_min_natural_value_, _t_max_natural_value_, _t_min_override_val_, _t_max_override_val_;
	bool _override_t_min_, _override_t_max_;
	mutable bool _record_full_data_;
	bool _host_loaded_, _satellite_loaded_;
	mutable bool _calculated_, _bad_result_, _current_satellite_in_use_,
			_current_host_in_use_;
	bool _evolving_host_;
	bool _using_private_init_host_, _using_private_init_satellite_;
	tNFW_profile _private_tNFW_init_host_, _private_tNFW_init_satellite_;

	// Splines for orbit data
	// Must be mutable since the spline class used doesn't allow conceptual constness for calculating
	mutable SALTSA::interpolator _x_spline_, _y_spline_, _z_spline_,
			_test_mass_spline_;
	mutable spline_derivative _vx_spline_, _vy_spline_, _vz_spline_;
	mutable std::vector< SALTSA::interpolator > _host_parameter_splines_;

	// Vectors for output data
	mutable std::vector< double > _delta_rho_list_, _sum_delta_rho_list_,
			_x_data_, _y_data_, _z_data_, _vx_data_, _vy_data_, _vz_data_,
			_rt_list_, _rt_ratio_list_;
	mutable std::vector< double > mret_list;
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
			const double &t_step, const bool silent = false ) const;
	const double _step_length_factor( const double & v, const double & r ) const;
	const double _rvir( const int index = 0 ) const;
#endif

public:

	// Constructors, destructors, and related operations
	stripping_orbit_segment();
	stripping_orbit_segment(
			const stripping_orbit_segment &other_orbit_spline );
	stripping_orbit_segment & operator=(
			const stripping_orbit_segment &other_orbit_spline );
	stripping_orbit_segment( const density_profile *host,
			const density_profile *satellite,
			const int init_resolution = 200 );
	virtual ~stripping_orbit_segment();
	stripping_orbit_segment *stripping_orbit_spline_clone() const;

	// Functions to add points to the splines
	const int add_point( const double &x, const double &y,
			const double &z, const double &t,
			double new_test_mass = 1 );
	const int add_point( const double &x, const double &y,
			const double &z, const double &vx,
			const double &vy, const double &vz, const double &t,
			double new_test_mass = 1 );
	const int add_x_point( const double &x, const double &t );
	const int add_y_point( const double &y, const double &t );
	const int add_z_point( const double &z, const double &t );
	const int add_vx_point( const double &vx, const double &t );
	const int add_vy_point( const double &vy, const double &t );
	const int add_vz_point( const double &vz, const double &t );
	const int add_unknown_vx_point( const double &t );
	const int add_unknown_vy_point( const double &t );
	const int add_unknown_vz_point( const double &t );
	const int add_test_mass_point( const double test_mass, const double &t );
	const int add_host_parameter_point( const unsigned int num_parameters,
			const std::vector< double > &parameters, const double &t,
			const bool silent = false );


	// Setting integration parameters
#if(1)
	const int set_resolution( const int new_spline_resolution,
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
	const int set_tidal_shocking_amplification( const double new_tidal_shocking_amplification,
			const bool silent=false );
	const int set_tidal_shocking_persistance( const double new_tidal_shocking_persistance,
			const bool silent=false );
	const int set_tidal_shocking_power( const double new_tidal_shocking_power,
			const bool silent=false );
#endif


	// Set initial/global parameters
	const int set_tNFW_init_satellite( const double &new_init_mvir0,
			const double z = 0, const double new_init_c = 0,
			const double new_init_tau = 0 );
	const int set_tNFW_host( const double &new_mvir0, const double z = 0,
			const double new_c = 0, const double new_tau = 0 );
	const int set_t_min( const double &new_t_min );
	const int set_t_max( const double &new_t_max );
	const int reset_t_min();
	const int reset_t_max();
	const int set_init_sum_deltarho( const double &new_init_sum_deltarho );
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
			true, const double mret_multiplier = 1,
			const bool silent = false ) const;

	// Get current state of object
	const unsigned int length() const;

	// Accessors
#if(1)

	// Integration parameters
#if(1)
	const int & spline_resolution() const {return _spline_resolution_;}
	const double & v_0() const {return _v_0_;}
	const double & r_0() const {return _r_0_;}
	const double & step_length_power() const {return _step_length_power_;}
	const double & step_factor_max() const {return _step_factor_max_;}
	const double & step_factor_min() const {return _step_factor_min_;}
#endif

	// Tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	const double & tidal_stripping_amplification() const {return _tidal_stripping_amplification_;}
	const double & tidal_stripping_deceleration() const {return _tidal_stripping_deceleration_;}
	const double & tidal_shocking_amplification() const {return _tidal_shocking_amplification_;}
	const double & tidal_shocking_persistance() const {return _tidal_shocking_persistance_;}
	const double & tidal_shocking_power() const {return _tidal_shocking_power_;}
#endif

	const bool & calculated() const {return _calculated_;};
	const bool & bad_result() const {return _bad_result_;};
	const density_profile * init_satellite_ptr() const {return _init_satellite_ptr_;};
	const density_profile * init_host_ptr() const {return _init_host_ptr_;};
	const double & t_min_natural_value() const {return _t_min_natural_value_;};
#endif

	// Calculation assistance functions
	const double tidal_strip_retained( const density_profile *host,
			const density_profile *satellite, const double &r,
			const double &vr, const double &vt,
			const double &time_step, const double &sum_rho = 0 ) const;
	const double get_rt( const density_profile *host,
			const density_profile *satellite, const double &r,
			const double &vr, const double &vt,
			const double &time_step, const double &sum_rho,
			const bool silent = false ) const;

	// Get final data (returns 1 on error)
	const int get_final_mret( double & mret,
			const bool silent = false ) const;
	const int get_final_sum_deltarho( double & final_sum_deltarho,
			const bool silent = false ) const;
	const int get_final_fmret( double & final_fmret,
			const bool silent = false ) const;
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt, const bool silent =
			false ) const;
	const int clone_final_satellite( density_profile * & final_satellite_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!
	const int clone_final_host( density_profile * & final_host_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!

	// Get final data (throws exception on error)
	const double final_mret() const;
	const double final_sum_deltarho() const;
	const double final_fmret() const;
	const gabdt final_sum_gabdt() const;
	const density_profile * final_satellite() const; // Creates a clone. Make sure to delete!
	const density_profile * final_host() const; // Creates a clone. Make sure to delete!

};
// class stripping_orbit_segment

class solve_rt_it_function: public functor< double > // Always uses one param, returns new value for that parameter for iteration.
{
	/************************************************************
	 solve_rt_it_function
	 --------------------

	 Child of function_class

	 Provides a function_class * to be used to iteratively solve
	 for tidal radius.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	double sum_rho, Daccel, omega;
	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;
	solve_rt_it_function( const double init_omega,
			const density_profile *init_satellite, const double init_Daccel,
			const double init_sum_rho = 0 );

	solve_rt_it_function();
};

class solve_rt_grid_function: public functor< double >
{
	/************************************************************
	 solve_rt_grid_function
	 ----------------------

	 Child of function_class

	 Provides a function_class * to be used with the grid solver.

	 \************************************************************/
public:
	const density_profile *satellite_ptr;

	double sum_rho, Daccel, omega;
	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;
	solve_rt_grid_function( const double init_omega,
			const density_profile *init_satellite, const double init_Daccel,
			const double init_sum_rho = 0 );

	solve_rt_grid_function();
};

#endif // end class definitions

} // end namespace SALTSA

#endif // __SALTSA_ORBIT_H_INCLUDED__
