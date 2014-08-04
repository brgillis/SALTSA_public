/**********************************************************************
 brg_orbit.h
 -----------

 This is the header file for classes and functions which handle tidal
 stripping along orbits.

 For the casual user, only the stripping_orbit class needs to be
 understood (and perhaps some fiddling with the Global Parameters in this
 file for tuning purposes). This class stores a number of points
 representing a satellite's position and velocity relative to its host
 halo at various times.

 In order to successfully calculate stripping, the class must be
 assigned a host density profile, an initial satellite density profile,
 and at least two points along the orbit. The density profiles should be
 derived classes of brgastro::density_profile, and pointers to them
 should be passed to this class using one of the method described below.

 When the calc() function is called (either by the user or by another
 function), the class performs the following operations:

 -Calculates splines representing the satellite's position and velocity
 over time, allowing these values to be estimated at points intermediate
 the known values.
 -Integrates along the path of the orbit. At each step, it:
 --Estimates the tidal radius and the mass loss
 --Determines how the satellite's density profile will change from this
 mass loss
 --Estimates the amount of energy injected due to tidal shocking, and
 converts this into an estimate of the effect decrease of the
 satellite's density.

 After this, output functions can be called without needing to
 recalculate, unless any of the input parameters are changed.


 Public methods for brgastro::stripping_orbit
 --------------------------------------------


 Constructors:


 brgastro::stripping_orbit::stripping_orbit()

 Creates an empty object.


 brgastro::stripping_orbit::stripping_orbit(
 const stripping_orbit &other_stripping_orbit)

 Creates as a copy. Pointers for host and initial satellite are copied
 if they point to external profiles, cloned if they point to private
 profiles.

 An assignment operator (=) is also defined, which performs largely the
 same function, and also handles cleanup of data being overwritten.


 brgastro::stripping_orbit::stripping_orbit(density_profile *host,
 density_profile *init_satellite, const int init_resolution = 200);

 Creates an object and assigns pointers to the density_profile for the
 host and initial satellite. It also optionally defines the resolution
 to be used for calculating stripping.


 Destructor:


 brgastro::stripping_orbit::~stripping_orbit()

 Performs cleanup on any dynamically-created objects. Declared as
 virtual in case it's overwritten in the future.


 Assignment:


 const int brgastro::stripping_orbit::set_init_satellite(
 density_profile *new_init_satellite)

 Sets a pointer to a halo which represents the initial satellite. The
 pointed-to object will not be altered by this class.


 const int brgastro::stripping_orbit::set_tNFW_init_satellite(
 double new_init_mvir0, double new_z=0, double new_init_c=0,
 double new_init_tau=0)

 Creates a private initial satellite to use, using the predefined
 tNFW_profile, from the passed parameters. mvir0 must be passed in
 units of kg.


 const int brgastro::stripping_orbit::set_host(density_profile *new_host)

 Sets a pointer to a halo which represents the host halo. The pointed-to
 object will not be altered by this class.


 const int brgastro::stripping_orbit::set_tNFW_host(double new_mvir0,
 double new_z=0, double new_c=0, double new_tau=0)

 Creates a private host profile to use, using the predefined
 tNFW_profile, from the passed parameters. mvir0 must be passed in
 units of kg.


 const int brgastro::stripping_orbit::add_point(
 double  x, double  y, double  z,
 double vx, double vy, double vz,
 unit_time t, double test_mass=1)

 Functions to add a point to the satellite's orbit. Distances must be
 given in units of m, velocities in m/s. Optionally, the user may also
 include a comparison "test mass" value for each point, representing the
 retained mass fraction at the timestep as determined by some other
 method. These values can later be printed alongside the values#
 calculated by this function for comparison.


 const int brgastro::stripping_orbit::clear()

 Clears all calculations and assigned values.


 const int brgastro::stripping_orbit::clear_calcs()

 Clears only calculated values, for instance if you want to calculate
 again with a different resolution.


 Calculation functions:


 const int brgastro::stripping_orbit::calc(const int resolution=0,
 const bool record_data=false)

 Calculates stripping for the orbit. If the resolution parameter is set
 here, it will override the default resolution used by the class. If
 record_data is set to true, the class will record data on the
 satellite's position, velocity, retained mass fraction, and halo
 parameters at each time step, allowing a data table to be output for
 inspection.


 const int brgastro::stripping_orbit::get_final_mret(
 BRG_MASS & mret )

 Writes the final retained mass into the passed variable. The type that
 must be passed depends on whether the _BRG_USE_UNITS_ compiler flag is
 defined or not (if it is defined, brgastro::unit_mass, else double).

 const BRG_MASS brgastro::stripping_orbit::final_mret()

 Similar to above, except returns the retained mass and throws an
 exception if the orbit's stripping cannot be calculated.


 const int brgastro::stripping_orbit::get_final_fmret(
 double & fmret )

 Writes the final fraction of mass retained into the passed variable.

 const double brgastro::stripping_orbit::final_fmret()

 Similar to above, except returns the fraction of retained mass and
 throws an exception if the orbit's stripping cannot be calculated.


 Detailed output functions:


 const int brgastro::stripping_orbit::set_satellite_output_parameters(
 const unsigned int num_parameters,
 std::vector< bool > satellite_output_parameters )

 A function to tell the class which of the satellite's density profile's
 parameters (from the get_parameters function) you wish to be output
 in the resulting data table. For instance, the tNFW_profile will give
 parameters of [ mvir0, z, c, tau ], but only tau will be altered by
 the simulated stripping. To output only tau, you would pass a vector of
 [false, false, false, true].


 const int brgastro::stripping_orbit::set_satellite_parameter_unitconvs(
 const unsigned int num_parameters, std::vector< double > satellite_unitconvs)

 Tells the class how to convert the units of the halo's output parameters
 before printing them. For instance, if you wanted a tNFW_profile's
 mvir0 in 10^10 Msun and its other parameters unchanged, you would pass
 [unitconv::ttMsuntokg, 1, 1, 1].


 const int brgastro::stripping_orbit::print_full_data( std::ostream *out)

 Prints a data table in ASCII format into the passed stream. This will
 only work if the stripping was calculated with the record_data parameter
 set to true. One row will be printed per time step (note that the
 step size is variable unless step_length_power is set to 0 in the
 configuration below).

 The default columns:

 "t"  - Time in Gyr
 "x"  - x position in kpc
 "y"  - y position in kpc
 "z"  - z position in kpc
 "vx" - x velocity in kpc/Gyr
 "vy" - y velocity in kpc/Gyr
 "vz" - z velocity in kpc/Gyr
 "d"  - radial distance from host in kpc
 "v"  - speed in kpc/Gyr
 "m_ret" - Retained mass fraction
 "m_frac_lost" - Rate of fractional mass loss (mass fraction loss over
 last time step divided by length of last time step)
 "sim_m_ret" - Retained mass fraction in comparison data
 "sim_m_frac_lost" - Rate of fractional mass loss in comparison data
 "rt" - Tidal radius of the satellite in kpc
 "rt/rvir" - Tidal radius of the satellite compared to its virial radius

 If the set_satellite_output_parameters function has been used,
 additional columns may be printed representing the parameters of the
 satellite halo at each time step.


 See the include brg_orbit_example.cpp file for an example on using this
 class.


 \**********************************************************************/

#ifndef __BRG_ORBIT_H_INCLUDED__
#define __BRG_ORBIT_H_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <stdexcept>
#include "brg_interpolator.h"
#include "brg_astro.h"
#include "brg_phase.hpp"
#include "brg_units.h"

namespace brgastro
{

/** Class Forward Declarations **/
#if (1)
// These declarations allow the various classes to point to each other without worry about
// which order they're declared in. (Order does still matter for classes containing other
// classes, though.)
//
// Classes are document further in the definitions
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

class interpolator_functor: public functor< BRG_UNITS >
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

	const brgastro::interpolator *_interpolator_ptr_;
	bool _interpolator_ptr_set_up_;

public:

	// Swap functions
	void swap(interpolator_functor & other);
	friend void swap(interpolator_functor & same, interpolator_functor & other) {same.swap(other);}

	// Constructors
	interpolator_functor();
	interpolator_functor(const interpolator_functor& other);
	interpolator_functor(const brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_functor()
	{
	}

	// Operator=
	interpolator_functor& operator=(interpolator_functor other);

	// Set functions
	const int set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	const int operator()( const BRG_UNITS & in_param,
			BRG_UNITS & out_param, const bool silent = false ) const;

};
// class interpolator_functor

class interpolator_derivative_functor: public functor< BRG_UNITS >
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
	interpolator_derivative_functor( brgastro::interpolator *init_interpolator_ptr );

	// Destructor
	virtual ~interpolator_derivative_functor()
	{
	}

	// Operator=
	interpolator_derivative_functor& operator=(interpolator_derivative_functor other);

	// Set functions
	const int set_interpolator_ptr( const brgastro::interpolator *new_interpolator_ptr );

	// Function method
	const int operator()( const BRG_UNITS & in_param,
			BRG_UNITS & out_param, const bool silent = false ) const;

};
// class interpolator_derivative_functor

class interpolator_derivative_weight_functor: public functor< BRG_UNITS >
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
	const int operator()( const BRG_UNITS & in_param,
			BRG_UNITS & out_param, const bool silent = false ) const;
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
	brgastro::interpolator *_interpolator_ptr_;
	mutable brgastro::interpolator _known_interpolator_, _estimated_interpolator_;
	bool _interpolator_ptr_set_up_;
	mutable bool _calculated_;

	interpolator_functor _interpolator_func_;

	std::vector< double > _unknown_t_list_;

	static double _default_sample_scale_, _default_sample_max_width_,
			_default_sample_precision_;
	double _sample_scale_, _sample_max_width_, _sample_precision_;
	mutable double _t_min_, _t_max_;

	brgastro::interpolator::allowed_interpolation_type _interpolation_type_;

public:

	// Swap functions
	void swap(interpolator_derivative &other);
	friend void swap(interpolator_derivative &same, interpolator_derivative &other) {same.swap(other);}

	// Constructors
	interpolator_derivative();
	interpolator_derivative( const interpolator_derivative &other);
	interpolator_derivative( brgastro::interpolator *init_interpolator_ptr );

	// Destructors
	virtual ~interpolator_derivative()
	{
	}

	// Operator=
	interpolator_derivative & operator=(interpolator_derivative other);

	// Set functions
	const int set_spline_ptr( brgastro::interpolator *new_interpolator_ptr );
	const int clear_spline_ptr();
	const int set_default_sample_scale( double new_default_sample_scale );
	const int set_default_sample_max_width(
			double new_default_sample_max_width );
	const int set_sample_scale( double new_sample_scale );
	const int set_sample_max_width( double new_sample_max_width );
	const int reset_sample_scale(); // Sets it to default
	const int reset_sample_max_width(); // Sets it to default
	const int set_interpolation_type(
			brgastro::interpolator::allowed_interpolation_type new_interpolation_type);
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

	BRG_DISTANCE _x_, _y_, _z_, _r_;
	BRG_TIME _dt_;
	mutable std::vector< std::vector< long double > > _dv_;

public:

	// Swap functions
	void swap(gabdt &other);
	friend void swap(gabdt &same, gabdt &other) {same.swap(other);}

	// Constructors
	gabdt();
	gabdt( const gabdt & other_gabdt );
	gabdt( const density_profile *init_host, const BRG_DISTANCE &init_x,
			const BRG_DISTANCE &init_y, const BRG_DISTANCE &init_z,
			const BRG_TIME &init_dt );

	// Destructor
	virtual ~gabdt();

	// Full clear function
	const int clear();

	// Set functions
	const int set( const brgastro::density_profile *new_host_ptr,
			const BRG_DISTANCE &new_x, const BRG_DISTANCE &new_y,
			const BRG_DISTANCE &new_z, const BRG_TIME &new_dt );
	const int set_pos( const BRG_DISTANCE &new_x, const BRG_DISTANCE &new_y,
			const BRG_DISTANCE &new_z );
	const int set_dt( const BRG_TIME &dt );
	const int set_host_ptr( const density_profile *new_host_ptr );
	const int override_zero();

	// Calculation function
	const int calc_dv( const bool silent = false ) const;

	// Get functions
	const density_profile * host() const;
	const BRG_DISTANCE x() const;
	const BRG_DISTANCE y() const;
	const BRG_DISTANCE z() const;
	const BRG_DISTANCE r() const;
	const std::vector< std::vector< long double > > dv() const; // involves calculation if necessary
	const long double dv( const int x_i, const int y_i ) const; // involves calculation if necessary

	// Operator overloading
	const BRG_UNITS operator*( const gabdt & other_gabdt ) const; // Dot-product(ish) operator

	gabdt & operator=( gabdt other_gabdt ); // Assignment

	gabdt & operator+=( const gabdt & other_gabdt ); // Addition
	gabdt operator+( const gabdt & other_gabdt ) const;

	gabdt & operator*=( const double scale_fraction ); // Multiplication by a double
	gabdt operator*( const double scale_fraction ) const;

};

class gabdt_functor: public functor< std::vector< BRG_UNITS > >
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
	const int operator()( const std::vector< BRG_UNITS > & in_params,
			std::vector< BRG_UNITS > & out_params,
			const bool silent = false ) const;

};

class stripping_orbit
{
	/************************************************************
	 stripping_orbit
	 --------------

	 This is the primary class used for calculating tidal
	 stripping. See documentation at the top of this file for a
	 full briefing on what you should do with it.

	 \************************************************************/
public:
	// Define the allowed interpolation types with an enum

	enum allowed_interpolation_type {
		LOWER,
		UPPER,
		LINEAR,
		SPLINE,
		UNSET
	}; // end enum allowed_interpolation_type

private:
#if (1)

	// Default integration parameters
#if(1)
	// Default number of steps for which stripping is calculated
	static int _default_spline_resolution_;

	// Default interpolation method
	static allowed_interpolation_type _default_interpolation_type_;

	// Variable step length tweaking: Time step length is proportional to (v_0/v)^(step_length_power)
	// This gives smaller steps when the satellite is moving faster.
	// If you want to turn off adaptive step size, set step_length_power to 0
	// Alternatively, set step_length_power to 1 for even steps in position
	static double _default_v_0_; // 400 km/s
	static double _default_r_0_; // 400 kpc
	static double _default_step_length_power_;
	static double _default_step_factor_max_; // Maximum allowed value of (v_0/v)^(step_length_power)
	static double _default_step_factor_min_; // Minimum allowed value of (v_0/v)^(step_length_power)
#endif

	// Integration parameters
#if(1)
	// Number of steps for which stripping is calculated
	int _spline_resolution_;

	// Interpolation method
	allowed_interpolation_type _interpolation_type_;

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

	// Default tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	static double _default_tidal_stripping_amplification_; // Amplifies tidal stripping by this factor
	static double _default_tidal_stripping_deceleration_; // If positive, increase tidal stripping near pericentre,
														  // if negative, decrease near pericentre
	static double _default_tidal_shocking_amplification_; // Amplifies tidal heating by this factor
	static double _default_tidal_shocking_persistance_; // How long shocking is active for
	static double _default_tidal_shocking_power_; // Affects interplay of stripping and satellite halo profile
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

	// Global info for the orbit
#if(1)
	mutable int _num_segments_;
	BRG_TIME _t_min_natural_value_, _t_max_natural_value_,
	         _t_min_override_value_, _t_max_override_value_;
	bool _override_t_min_, _override_t_max_; // Tells if min/max have been set manually, so the manual settings can be used
	mutable bool _likely_disrupted_;
#endif

	// Data for output info for satellite and host_ptr
#if(1)
	std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;
#endif

	// Lists of points on the orbit and related info
#if(1)
	std::vector< std::pair< double, double > > _x_points_,
			_y_points_, _z_points_, _test_mass_points_;
	std::vector< std::pair< double, double > > _vx_points_,
			_vy_points_, _vz_points_;
	std::vector< double > _vx_unknown_points_,
			_vy_unknown_points_, _vz_unknown_points_,
			_t_points_, _host_param_t_points_;
	std::vector< std::pair< double, std::vector< BRG_UNITS > > > _host_parameter_points_;
	std::vector< double > _discontinuity_times_;
	mutable std::vector< double > _cleaned_discontinuity_times_;
	mutable brgastro::interpolator _m_ret_interpolator_;
	mutable brgastro::interpolator _test_mass_interpolator_, _test_mass_error_interpolator_;
	mutable std::vector< double > _final_fmret_list_;
	mutable std::vector< brgastro::stripping_orbit_segment > _orbit_segments_;
#endif

	// Host and satellite pointers and info
#if(1)
	const density_profile *_init_host_ptr_, *_init_satellite_ptr_;
	mutable bool _record_full_data_;
	bool _host_is_evolving_;
	bool _host_loaded_, _satellite_loaded_;
	mutable bool _calculated_, _bad_result_;
	bool _using_private_init_host_, _using_private_init_satellite_;
	tNFW_profile _private_tNFW_init_host_, _private_tNFW_init_satellite_;
	mutable std::vector< brgastro::stripping_orbit_segment >::iterator _final_good_segment_;
#endif

	const std::vector< brgastro::stripping_orbit_segment >::iterator _final_good_segment() const;

	// Private methods
	const int _pass_parameters_to_segment(
			brgastro::stripping_orbit_segment & segment,
			brgastro::density_profile *temp_satellite=NULL,
			brgastro::density_profile *temp_host=NULL,
			unsigned int resolution=0) const;
#endif
public:
#if (1)

	// Swap functions
	void swap(stripping_orbit &other);
	friend void swap(stripping_orbit &same, stripping_orbit &other) {same.swap(other);}

	// Constructors, destructor, and related functions
	stripping_orbit(); // Basic constructor
	stripping_orbit( const density_profile *host, const density_profile *satellite,
			const int init_resolution = 200 ); // Constructor with initial params
	stripping_orbit( const stripping_orbit &other_orbit_spline ); // Copy constructor
	stripping_orbit & operator=( stripping_orbit other ); // Assignment operator
	stripping_orbit *stripping_orbit_clone(); // Clone function
	virtual ~stripping_orbit(); // Virtual destructor

	// Setting default integration parameters
#if(1)
	static const int set_default_resolution( const int new_default_spline_resolution);
	const int set_default_resolution( const int new_default_spline_resolution,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_interpolation_type(
			const allowed_interpolation_type new_default_interpolation_type);
	const int set_default_interpolation_type(
			const allowed_interpolation_type new_default_interpolation_type,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_v_0( const double new_default_v_0);
	const int set_default_v_0( const double new_default_v_0,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_r_0( const double new_default_r_0);
	const int set_default_r_0( const double new_default_r_0,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_step_length_power( const double new_default_step_length_power);
	const int set_default_step_length_power( const double new_default_step_length_power,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_step_factor_max( const double new_default_step_factor_max);
	const int set_default_step_factor_max( const double new_default_step_factor_max,
			const bool override_current=false,
			const bool silent=false );
	static const int set_default_step_factor_min( const double new_default_step_factor_min);
	const int set_default_step_factor_min( const double new_default_step_factor_min,
			const bool override_current=false,
			const bool silent=false );
#endif

	// Setting default tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	static const int set_default_tidal_stripping_amplification(
			const double new_default_tidal_stripping_amplification);
	const int set_default_tidal_stripping_amplification(
			const double new_default_tidal_stripping_amplification,
			const bool override_current,
			const bool silent=false );
	static const int set_default_tidal_stripping_deceleration(
			const double new_default_tidal_stripping_deceleration);
	const int set_default_tidal_stripping_deceleration(
			const double new_default_tidal_stripping_deceleration,
			const bool override_current,
			const bool silent=false );
	static const int set_default_tidal_shocking_amplification(
			const double new_default_tidal_shocking_amplification);
	const int set_default_tidal_shocking_amplification(
			const double new_default_tidal_shocking_amplification,
			const bool override_current,
			const bool silent=false );
	static const int set_default_tidal_shocking_persistance(
			const double new_default_tidal_shocking_persistance);
	const int set_default_tidal_shocking_persistance(
			const double new_default_tidal_shocking_persistance,
			const bool override_current,
			const bool silent=false );
	static const int set_default_tidal_shocking_power(
			const double new_default_tidal_shocking_power);
	const int set_default_tidal_shocking_power(
			const double new_default_tidal_shocking_power,
			const bool override_current,
			const bool silent=false );
#endif

	// Setting integration parameters
#if(1)
	const int set_resolution( const int new_spline_resolution,
			const bool silent=false );
	const int set_interpolation_type( const allowed_interpolation_type new_type,
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

	// Resetting integration parameters
#if(1)
	const int reset_resolution();
	const int reset_interpolation_type();
	const int reset_v_0();
	const int reset_r_0();
	const int reset_step_length_power();
	const int reset_step_factor_max();
	const int reset_step_factor_min();
#endif

	// Resetting tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	const int reset_tidal_stripping_amplification();
	const int reset_tidal_stripping_deceleration();
	const int reset_tidal_shocking_amplification();
	const int reset_tidal_shocking_persistance();
	const int reset_tidal_shocking_power();
#endif

	// Adding data points and clearing those vectors
	const int add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_VELOCITY &vx,
			const BRG_VELOCITY &vy, const BRG_VELOCITY &vz, const BRG_TIME &t,
			const double test_mass = 1, const double test_mass_error = 1 );
	const int force_add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_VELOCITY &vx,
			const BRG_VELOCITY &vy, const BRG_VELOCITY &vz, const BRG_TIME &t,
			const double test_mass = 1, const double test_mass_error = 1 );
	const int add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_TIME &t,
			const double new_test_mass = 1, const double test_mass_error = 1 ); // Only use if v is unknown
	const int force_add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_TIME &t,
			const double new_test_mass = 1, const double test_mass_error = 1 ); // Only use if v is unknown
	const int add_discontinuity_time( const BRG_TIME &t ); // Splits into segments to be calculated individually
	const int add_host_parameter_point( const std::vector< BRG_UNITS > &parameters, const BRG_TIME &t,
			const bool silent = false ); // Tells how host_ptr is evolving
	const int force_add_host_parameter_point( const std::vector< BRG_UNITS > &parameters, const BRG_TIME &t,
			const bool silent = false ); // Tells how host_ptr is evolving
	const int clear_points();
	const int clear_discontinuity_times();
	const int clear_host_parameter_points();

	// Setting and clearing _init satellite/host_ptr
	const int set_init_satellite( const density_profile *new_init_satellite ); // Uses pointer to existing profile
	const int set_init_host( const density_profile *new_init_host ); // Uses pointer to existing profile
	const int set_tNFW_init_satellite( const BRG_MASS &new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 ); // Creates a new profile
	const int set_tNFW_init_host( const BRG_MASS &new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 ); // Creates a new profile
	const int clear_init_satellite();
	const int clear_init_host();

	// Setting and resetting _t_min_/max
	const int set_t_min( const BRG_TIME &new_t_min ); // In case you want to integrate only part of the submitted orbit, or go outside the bounds
	const int set_t_max( const BRG_TIME &new_t_max ); // In case you want to integrate only part of the submitted orbit, or go outside the bounds
	const int reset_t_min();
	const int reset_t_max();

	// Functions for determining how calc() will be called
	const int set_record_full_data( const bool new_record_full_data ) const;

	// Global clearing functions
	const int clear();
	const int clear_calcs() const; // Only clears calculated (mutable) values

	// Function to force calculation
	const int calc( const bool silent = false ) const; // Using conceptual const-ness

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

	// Print output functions
	const int print_full_data( std::ostream *out ) const;
	const int print_segment_data( std::ostream *out,
			const int segment_number ) const;

	// Accessors to private data
#if(1)

	// Default integration parameters
#if(1)
	static const int & default_spline_resolution() {return _default_spline_resolution_;}
	static const allowed_interpolation_type & default_interpolation_type()
		{return _default_interpolation_type_;}
	static const double & default_v_0() {return _default_v_0_;}
	static const double & default_r_0() {return _default_r_0_;}
	static const double & default_step_length_power() {return _default_step_length_power_;}
	static const double & default_step_factor_max() {return _default_step_factor_max_;}
	static const double & default_step_factor_min() {return _default_step_factor_min_;}
#endif

	// Default tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	static const double & default_tidal_stripping_amplification() {return _default_tidal_stripping_amplification_;}
	static const double & default_tidal_stripping_deceleration() {return _default_tidal_stripping_deceleration_;}
	static const double & default_tidal_shocking_amplification() {return _default_tidal_shocking_amplification_;}
	static const double & default_tidal_shocking_persistance() {return _default_tidal_shocking_persistance_;}
	static const double & default_tidal_shocking_power() {return _default_tidal_shocking_power_;}
#endif


	// Integration parameters
#if(1)
	const int & spline_resolution() const {return _spline_resolution_;}
	const allowed_interpolation_type & interpolation_type()
		{return _interpolation_type_;}
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

	const int & num_segments() const {return _num_segments_;};
	const BRG_TIME & t_min_natural_value() const {return _t_min_natural_value_;};
	const BRG_TIME & t_max_natural_value() const {return _t_max_natural_value_;};
	const BRG_TIME & t_min_override_value() const {return _t_min_override_value_;};
	const BRG_TIME & t_max_override_value() const {return _t_max_override_value_;};
	const BRG_TIME & t_min() const
	{
		if(_override_t_min_)
			return _t_min_override_value_;
		else
			return _t_min_natural_value_;
	}
	const BRG_TIME & t_max() const
	{
		if(_override_t_max_)
			return _t_max_override_value_;
		else
			return _t_max_natural_value_;
	}
	const bool & override_t_min() const {return _override_t_min_;};
	const bool & override_t_max() const {return _override_t_max_;};

	const std::vector< double > & satellite_parameter_unitconvs() const {return _satellite_parameter_unitconvs_;};
	const std::vector< double > & host_parameter_unitconvs() const {return _host_parameter_unitconvs_;};
	const std::vector< bool > & satellite_output_parameters() const {return _satellite_output_parameters_;};
	const std::vector< bool > & host_output_parameters() const {return _host_output_parameters_;};

	const std::vector< std::pair< double, double > > & x_spline_points() const {return _x_points_;};
	const std::vector< std::pair< double, double > > & y_spline_points() const {return _y_points_;};
	const std::vector< std::pair< double, double > > & z_spline_points() const {return _z_points_;};
	const std::vector< std::pair< double, double > > & vx_spline_points() const {return _vx_points_;};
	const std::vector< std::pair< double, double > > & vy_spline_points() const {return _vy_points_;};
	const std::vector< std::pair< double, double > > & vz_spline_points() const {return _vz_points_;};
	const std::vector< double > & vx_spline_unknown_points() const {return _vx_unknown_points_;};
	const std::vector< double > & vy_spline_unknown_points() const {return _vy_unknown_points_;};
	const std::vector< double > & vz_spline_unknown_points() const {return _vz_unknown_points_;};
	const std::vector< std::pair< double, double > > & test_mass_spline_points() const {return _test_mass_points_;};
	const std::vector< std::pair< double, std::vector< BRG_UNITS > > > & host_parameter_spline_points() const {return _host_parameter_points_;};
	const std::vector< double > & discontinuity_times() const {return _discontinuity_times_;};

	const density_profile * init_satellite_ptr() const {return _init_satellite_ptr_;};
	const density_profile * init_host_ptr() const {return _init_host_ptr_;};

	const bool & record_full_data() const {return _record_full_data_;};
	const bool & host_is_evolving() const {return _host_is_evolving_;};
	const bool & host_loaded() const {return _host_loaded_;};
	const bool & satellite_loaded() const {return _satellite_loaded_;};
	const bool & calculated() const {return _calculated_;};
	const bool & bad_result() const {return _bad_result_;};
	const bool & using_private_init_host() const {return _using_private_init_host_;};
	const bool & using_private_init_satellite() const {return _using_private_init_satellite_;};

	const tNFW_profile & private_tNFW_init_host() const {return _private_tNFW_init_host_;};
	const tNFW_profile & private_tNFW_init_satellite() const {return _private_tNFW_init_satellite_;};

	const std::vector<stripping_orbit_segment> & orbit_segments() const {return _orbit_segments_;};
#endif

	// Getting resultant data
#if (1)

	// Get final data (returns 1 on failure)
	const int get_final_mret( BRG_MASS & mret ) const;
	const int get_final_sum_deltarho( long double & final_sum_deltarho ) const;
	const int get_final_sum_deltarho( BRG_UNITS & final_sum_deltarho ) const;
	const int get_final_fmret( double & final_fmret ) const;
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt ) const;
	const int get_last_infall_time( BRG_TIME & t ) const;
	const int clone_final_satellite(
			density_profile * & final_satellite_clone ) const; // Creates clone. Make sure to delete!
	const int clone_final_host( density_profile * & final_host_clone ) const; // Creates clone. Make sure to delete!

	// Get final data (throws exception on failure)
	const BRG_MASS final_mret() const;
	const BRG_UNITS final_sum_deltarho() const;
	const double final_fmret() const;
	const gabdt final_sum_gabdt() const;
	const BRG_TIME last_infall_time() const;
	const bool & likely_disrupted() const;
	const density_profile * final_satellite() const; // Creates clone. Make sure to delete!
	const density_profile * final_host() const; // Creates clone. Make sure to delete!

	// Get data at arbitrary time (returns 1 on failure)
	const int get_mret_at_t( const BRG_TIME & t, BRG_MASS & mret) const; // Requires _record_full_data_
	const int get_fmret_at_t( const BRG_TIME & t, double & fmret ) const; // Requires _record_full_data_
	const int get_comp_fmret_at_t( const BRG_TIME & t, double & mret) const;
	const int get_comp_fmret_error_at_t( const BRG_TIME & t, double & mret) const;

	// Get data at arbitrary time (throws exception on failure)
	const BRG_MASS mret_at_t(const BRG_TIME & t) const; // Requires _record_full_data_
	const double fmret_at_t(const BRG_TIME & t) const; // Requires _record_full_data_
	const double comp_fmret_at_t(const BRG_TIME & t) const;
	const double comp_fmret_error_at_t(const BRG_TIME & t) const;

	// Quality of fit. Both versions require _record_full_data_
	const int get_quality_of_fit( double & Q, const unsigned int samples=100);
	const double quality_of_fit(const unsigned int samples=100);

#endif // Getting resultant data

#endif

};
// class stripping_orbit

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
	mutable std::vector< long double > _sum_delta_rho_list_, _mret_list_;
	mutable std::vector< std::vector< BRG_UNITS > > _satellite_parameter_data_; // Keeps track of satellite's parameters (ie. mass, tau)
	mutable std::vector< std::vector< BRG_UNITS > > _host_parameter_data_;

	// Other maintained vectors for calculations
	mutable std::vector< gabdt > _gabdt_list_;
	mutable std::vector< gabdt > _sum_gabdt_list_;
	mutable std::vector< brgastro::phase > _phase_list_, _phase_output_list_;

	// Output-modifying parameters
	int _num_parameters_;
	mutable std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;

	// Private functions
	const int _init();
	const int _reserve( int n, const bool silent = false ) const;
	const BRG_UNITS _delta_rho( const int index, const double x,
			const BRG_TIME &t_step, const bool silent = false ) const;
	const double _step_length_factor( const BRG_VELOCITY & v, const BRG_DISTANCE & r ) const;
	const BRG_DISTANCE _rvir( const int index = 0 ) const;
	const int _pass_interpolation_type() const;

	// Calculation assistance functions
	const double _tidal_strip_retained( const density_profile *host,
			const density_profile *satellite, const BRG_DISTANCE &r,
			const BRG_VELOCITY &vr, const BRG_VELOCITY &vt,
			const BRG_TIME &time_step, const long double &sum_delta_rho = 0 ) const;
	const BRG_DISTANCE _get_rt( const density_profile *host,
			const density_profile *satellite, const BRG_DISTANCE &r,
			const BRG_VELOCITY &vr, const BRG_VELOCITY &vt,
			const BRG_TIME &time_step, const long double &sum_delta_rho,
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
	const int add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_TIME &t,
			double new_test_mass = 1 );
	const int add_point( const BRG_DISTANCE &x, const BRG_DISTANCE &y,
			const BRG_DISTANCE &z, const BRG_VELOCITY &vx,
			const BRG_VELOCITY &vy, const BRG_VELOCITY &vz, const BRG_TIME &t,
			double new_test_mass = 1 );
	const int add_x_point( const BRG_DISTANCE &x, const BRG_TIME &t );
	const int add_y_point( const BRG_DISTANCE &y, const BRG_TIME &t );
	const int add_z_point( const BRG_DISTANCE &z, const BRG_TIME &t );
	const int add_vx_point( const BRG_VELOCITY &vx, const BRG_TIME &t );
	const int add_vy_point( const BRG_VELOCITY &vy, const BRG_TIME &t );
	const int add_vz_point( const BRG_VELOCITY &vz, const BRG_TIME &t );
	const int add_unknown_vx_point( const BRG_TIME &t );
	const int add_unknown_vy_point( const BRG_TIME &t );
	const int add_unknown_vz_point( const BRG_TIME &t );
	const int add_test_mass_point( const double test_mass, const BRG_TIME &t );
	const int add_host_parameter_point( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > &parameters, const BRG_TIME &t,
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
	const int set_tidal_shocking_amplification( const double new_tidal_shocking_amplification,
			const bool silent=false );
	const int set_tidal_shocking_persistance( const double new_tidal_shocking_persistance,
			const bool silent=false );
	const int set_tidal_shocking_power( const double new_tidal_shocking_power,
			const bool silent=false );
#endif


	// Set initial/global parameters
	const int set_tNFW_init_satellite( const BRG_MASS &new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 );
	const int set_tNFW_host( const BRG_MASS &new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 );
	const int set_t_min( const BRG_TIME &new_t_min );
	const int set_t_max( const BRG_TIME &new_t_max );
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
			true, const double mret_multiplier = 1,
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
	const BRG_TIME & t_min_natural_value() const {return _t_min_natural_value_;};
#endif

	// Get final data (returns 1 on error)
	const int get_final_mret( BRG_MASS & mret,
			const bool silent = false ) const;
	const int get_final_sum_deltarho( long double & final_sum_deltarho,
			const bool silent = false ) const;
	const int get_final_sum_deltarho( BRG_UNITS & final_sum_deltarho,
			const bool silent = false ) const;
	const int get_final_fmret( double & final_fmret,
			const bool silent = false ) const;
	const int get_mret_points( std::vector< std::pair<double,double> > & mret_points,
			const bool silent = false ) const;
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt, const bool silent =
			false ) const;
	const int clone_final_satellite( density_profile * & final_satellite_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!
	const int clone_final_host( density_profile * & final_host_clone,
			const bool silent = false ) const; // Creates a clone. Make sure to delete!

	// Get final data (throws exception on error)
	const BRG_MASS final_mret() const;
	const long double final_sum_deltarho() const;
	const double final_fmret() const;
	const std::vector< std::pair<double,double> > mret_points() const;
	const gabdt final_sum_gabdt() const;
	const bool & likely_disrupted() const;
	const density_profile * final_satellite() const; // Creates a clone. Make sure to delete!
	const density_profile * final_host() const; // Creates a clone. Make sure to delete!

};
// class stripping_orbit_segment

class solve_rt_it_functor: public functor< BRG_UNITS > // Always uses one param, returns new value for that parameter for iteration.
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
	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;
	solve_rt_it_functor( const BRG_UNITS init_omega,
			const density_profile *init_satellite, const BRG_UNITS init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_it_functor();
};

class solve_rt_grid_functor: public functor< BRG_UNITS >
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
	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;
	solve_rt_grid_functor( const BRG_UNITS init_omega,
			const density_profile *init_satellite, const BRG_UNITS init_Daccel,
			const long double init_sum_delta_rho = 0 );

	solve_rt_grid_functor();
};

#endif // end class definitions

} // end namespace brgastro

#endif //__BRG_ORBIT_H_INCLUDED__

