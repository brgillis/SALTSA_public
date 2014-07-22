/**********************************************************************
 @file SALTSA.h
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
 derived classes of SALTSA::density_profile, and pointers to them
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


 Public methods for SALTSA::stripping_orbit
 --------------------------------------------


 Constructors:


 SALTSA::stripping_orbit::stripping_orbit()

 Creates an empty object.


 SALTSA::stripping_orbit::stripping_orbit(
 const stripping_orbit &other_stripping_orbit)

 Creates as a copy. Pointers for host and initial satellite are copied
 if they point to external profiles, cloned if they point to private
 profiles.

 An assignment operator (=) is also defined, which performs largely the
 same function, and also handles cleanup of data being overwritten.


 SALTSA::stripping_orbit::stripping_orbit(density_profile *host,
 density_profile *init_satellite, const int init_resolution = 200);

 Creates an object and assigns pointers to the density_profile for the
 host and initial satellite. It also optionally defines the resolution
 to be used for calculating stripping.


 Destructor:


 SALTSA::stripping_orbit::~stripping_orbit()

 Performs cleanup on any dynamically-created objects. Declared as
 virtual in case it's overwritten in the future.


 Assignment:


 const int SALTSA::stripping_orbit::set_init_satellite(
 density_profile *new_init_satellite)

 Sets a pointer to a halo which represents the initial satellite. The
 pointed-to object will not be altered by this class.


 const int SALTSA::stripping_orbit::set_tNFW_init_satellite(
 double new_init_mvir0, double new_z=0, double new_init_c=0,
 double new_init_tau=0)

 Creates a private initial satellite to use, using the predefined
 tNFW_profile, from the passed parameters. mvir0 must be passed in
 units of kg.


 const int SALTSA::stripping_orbit::set_host(density_profile *new_host)

 Sets a pointer to a halo which represents the host halo. The pointed-to
 object will not be altered by this class.


 const int SALTSA::stripping_orbit::set_tNFW_host(double new_mvir0,
 double new_z=0, double new_c=0, double new_tau=0)

 Creates a private host profile to use, using the predefined
 tNFW_profile, from the passed parameters. mvir0 must be passed in
 units of kg.


 const int SALTSA::stripping_orbit::add_point(
 double  x, double  y, double  z,
 double vx, double vy, double vz,
 unit_time t, double test_mass=1)

 Functions to add a point to the satellite's orbit. Distances must be
 given in units of m, velocities in m/s. Optionally, the user may also
 include a comparison "test mass" value for each point, representing the
 retained mass fraction at the timestep as determined by some other
 method. These values can later be printed alongside the values#
 calculated by this function for comparison.


 const int SALTSA::stripping_orbit::clear()

 Clears all calculations and assigned values.


 const int SALTSA::stripping_orbit::clear_calcs()

 Clears only calculated values, for instance if you want to calculate
 again with a different resolution.


 Calculation functions:


 const int SALTSA::stripping_orbit::calc(const int resolution=0,
 const bool record_data=false)

 Calculates stripping for the orbit. If the resolution parameter is set
 here, it will override the default resolution used by the class. If
 record_data is set to true, the class will record data on the
 satellite's position, velocity, retained mass fraction, and halo
 parameters at each time step, allowing a data table to be output for
 inspection.


 const int SALTSA::stripping_orbit::get_final_mret(
 double & mret )

 Writes the final retained mass into the passed variable.

 const double SALTSA::stripping_orbit::final_mret()

 Similar to above, except returns the retained mass and throws an
 exception if the orbit's stripping cannot be calculated.


 const int SALTSA::stripping_orbit::get_final_fmret(
 double & fmret )

 Writes the final fraction of mass retained into the passed variable.

 const double SALTSA::stripping_orbit::final_fmret()

 Similar to above, except returns the fraction of retained mass and
 throws an exception if the orbit's stripping cannot be calculated.


 Detailed output functions:


 const int SALTSA::stripping_orbit::set_satellite_output_parameters(
 const unsigned int num_parameters,
 std::vector< bool > satellite_output_parameters )

 A function to tell the class which of the satellite's density profile's
 parameters (from the get_parameters function) you wish to be output
 in the resulting data table. For instance, the tNFW_profile will give
 parameters of [ mvir0, z, c, tau ], but only tau will be altered by
 the simulated stripping. To output only tau, you would pass a vector of
 [false, false, false, true].


 const int SALTSA::stripping_orbit::set_satellite_parameter_unitconvs(
 const unsigned int num_parameters, std::vector< double > satellite_unitconvs)

 Tells the class how to convert the units of the halo's output parameters
 before printing them. For instance, if you wanted a tNFW_profile's
 mvir0 in 10^10 Msun and its other parameters unchanged, you would pass
 [unitconv::ttMsuntokg, 1, 1, 1].


 const int SALTSA::stripping_orbit::print_full_data( std::ostream *out)

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

#ifndef __SALTSA_H_INCLUDED__
#define __SALTSA_H_INCLUDED__

#include <vector>
#include <cmath>
#include <iostream>

#include "SALTSA_global.h"

#include "SALTSA_interpolator.h"
#include "SALTSA_astro.h"

namespace SALTSA
{

 /** Class Forward Declarations **/
 #if (1)
 // These declarations allow the various classes to point to each other without worry about
 // which order they're declared in. (Order does still matter for classes containing other
 // classes, though.)
 class stripping_orbit_segment;
 class gabdt;

 #endif // end Class forward declarations

/** Class Definitions **/
#if (1)


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
	double _t_min_natural_value_, _t_max_natural_value_,
	         _t_min_override_value_, _t_max_override_value_;
	bool _override_t_min_, _override_t_max_; // Tells if min/max have been set manually, so the manual settings can be used
	mutable bool _likely_disrupted_;
#endif

	// Data for output info for satellite and host_ptr
#if(1)
	int _num_patameters_;
	std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;
#endif

	// Lists of points on the orbit and related info
#if(1)
	std::vector< std::pair< double, double > > _x_spline_points_,
			_y_spline_points_, _z_spline_points_, _d_spline_points_,
			_test_mass_spline_points_;
	std::vector< std::pair< double, double > > _vx_spline_points_,
			_vy_spline_points_, _vz_spline_points_;
	std::vector< double > _vx_spline_unknown_points_,
			_vy_spline_unknown_points_, _vz_spline_unknown_points_;
	std::vector< std::pair< double, std::vector< double > > > _host_parameter_spline_points_;
	std::vector< double > _discontinuity_times_;
	mutable std::vector< double > _cleaned_discontinuity_times_;
	mutable std::vector< double > _final_fmret_list_;
	int _num_discontinuities_;
	mutable int _num_cleaned_discontinuities_;
	mutable std::vector< SALTSA::stripping_orbit_segment > _orbit_segments_;
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
	mutable std::vector< SALTSA::stripping_orbit_segment >::iterator _final_good_segment_;
#endif

	const std::vector< SALTSA::stripping_orbit_segment >::iterator _final_good_segment() const;

	// Private methods
	const int _init();
	const int _pass_parameters_to_segment(
			SALTSA::stripping_orbit_segment & segment,
			SALTSA::density_profile *temp_satellite=NULL,
			SALTSA::density_profile *temp_host=NULL,
			unsigned int resolution=0) const;
#endif
public:
#if (1)
	// Constructors, destructor, and related functions
	stripping_orbit(); // Basic constructor
	stripping_orbit( const density_profile *host, const density_profile *satellite,
			const int init_resolution = 200 ); // Constructor with initial params
	stripping_orbit( const stripping_orbit &other_orbit_spline ); // Copy constructor
	stripping_orbit & operator=( const stripping_orbit &other_orbit_spline ); // Assignment operator
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
	const int set_default_tidal_stripping_amplification(
			const double new_default_tidal_stripping_amplification,
			const bool override_current=false,
			const bool silent=false );
	const int set_default_tidal_stripping_deceleration(
			const double new_default_tidal_stripping_deceleration,
			const bool override_current=false,
			const bool silent=false );
	const int set_default_tidal_shocking_amplification(
			const double new_default_tidal_shocking_amplification,
			const bool override_current=false,
			const bool silent=false );
	const int set_default_tidal_shocking_persistance(
			const double new_default_tidal_shocking_persistance,
			const bool override_current=false,
			const bool silent=false );
	const int set_default_tidal_shocking_power(
			const double new_default_tidal_shocking_power,
			const bool override_current=false,
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
	const int add_point( const double &x, const double &y,
			const double &z, const double &vx,
			const double &vy, const double &vz, const double &t,
			const double new_test_mass = 1 );
	const int add_point( const double &x, const double &y,
			const double &z, const double &t,
			const double new_test_mass = 1 ); // Only use if v is unknown
	const int add_discontinuity_time( const double &t ); // Splits into segments to be calculated individually
	const int add_host_parameter_point( const std::vector< double > &parameters, const double &t,
			const bool silent = false ); // Tells how host_ptr is evolving
	const int clear_points();
	const int clear_discontinuity_times();
	const int clear_host_parameter_points();

	// Setting and clearing _init satellite/host_ptr
	const int set_init_satellite( const density_profile *new_init_satellite ); // Uses pointer to existing profile
	const int set_init_host( const density_profile *new_init_host ); // Uses pointer to existing profile
	const int set_tNFW_init_satellite( const double &new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 ); // Creates a new profile
	const int set_tNFW_init_host( const double &new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 ); // Creates a new profile
	const int clear_init_satellite();
	const int clear_init_host();

	// Setting and resetting _t_min_/max
	const int set_t_min( const double &new_t_min ); // In case you want to integrate only part of the submitted orbit, or go outside the bounds
	const int set_t_max( const double &new_t_max ); // In case you want to integrate only part of the submitted orbit, or go outside the bounds
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
	const double & t_min_natural_value() const {return _t_min_natural_value_;};
	const double & t_max_natural_value() const {return _t_max_natural_value_;};
	const double & t_min_override_value() const {return _t_min_override_value_;};
	const double & t_max_override_value() const {return _t_max_override_value_;};
	const bool & override_t_min() const {return _override_t_min_;};
	const bool & override_t_max() const {return _override_t_max_;};

	const std::vector< double > & satellite_parameter_unitconvs() const {return _satellite_parameter_unitconvs_;};
	const std::vector< double > & host_parameter_unitconvs() const {return _host_parameter_unitconvs_;};
	const std::vector< bool > & satellite_output_parameters() const {return _satellite_output_parameters_;};
	const std::vector< bool > & host_output_parameters() const {return _host_output_parameters_;};

	const std::vector< std::pair< double, double > > & x_spline_points() const {return _x_spline_points_;};
	const std::vector< std::pair< double, double > > & y_spline_points() const {return _y_spline_points_;};
	const std::vector< std::pair< double, double > > & z_spline_points() const {return _z_spline_points_;};
	const std::vector< std::pair< double, double > > & vx_spline_points() const {return _vx_spline_points_;};
	const std::vector< std::pair< double, double > > & vy_spline_points() const {return _vy_spline_points_;};
	const std::vector< std::pair< double, double > > & vz_spline_points() const {return _vz_spline_points_;};
	const std::vector< double > & vx_spline_unknown_points() const {return _vx_spline_unknown_points_;};
	const std::vector< double > & vy_spline_unknown_points() const {return _vy_spline_unknown_points_;};
	const std::vector< double > & vz_spline_unknown_points() const {return _vz_spline_unknown_points_;};
	const std::vector< std::pair< double, double > > & d_spline_points() const {return _d_spline_points_;};
	const std::vector< std::pair< double, double > > & test_mass_spline_points() const {return _test_mass_spline_points_;};
	const std::vector< std::pair< double, std::vector< double > > > & host_parameter_spline_points() const {return _host_parameter_spline_points_;};
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

	// Get final data (returns 1 on failure)
	const int get_final_mret( double & mret ) const;
	const int get_final_sum_deltarho( long double & final_sum_deltarho ) const;
	const int get_final_sum_deltarho( double & final_sum_deltarho ) const;
	const int get_final_fmret( double & final_fmret ) const;
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt ) const;
	const int get_last_infall_time( double & t ) const;
	const int clone_final_satellite(
			density_profile * & final_satellite_clone ) const; // Creates clone. Make sure to delete!
	const int clone_final_host( density_profile * & final_host_clone ) const; // Creates clone. Make sure to delete!

	// Get final data (throws exception on failure)
	const double final_mret() const;
	const double final_sum_deltarho() const;
	const double final_fmret() const;
	const gabdt final_sum_gabdt() const;
	const double last_infall_time() const;
	const bool & likely_disrupted() const;
	const density_profile * final_satellite() const; // Creates clone. Make sure to delete!
	const density_profile * final_host() const; // Creates clone. Make sure to delete!
#endif

};
// class stripping_orbit

#endif // end class definitions

} // end namespace SALTSA

#include "SALTSA_orbit.h"

#endif // __SALTSA_H_INCLUDED__
