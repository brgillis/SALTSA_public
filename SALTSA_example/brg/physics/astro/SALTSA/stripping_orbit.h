/**********************************************************************\
  @file stripping_orbit.h

 This is the primary header file for classes and functions which handle
 tidal stripping along orbits.

 See the documentation of the stripping_orbit class below and the
 documentation of its member functions in the public section of this
 class for use notes, and see the included SALTSA_example.cpp file for
 an example on using this class.

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

// body file: brg/physics/astro/SALTSA/stripping_orbit.cpp

#ifndef _STRIPPING_ORBIT_H_INCLUDED_
#define _STRIPPING_ORBIT_H_INCLUDED_

#include <cstdlib>
#include <cmath>
#include <vector>
#include <stdexcept>

#include "brg/global.h"

#include "brg/math/interpolator/interpolator.h"
#include "brg/math/interpolator/interpolator_derivative.h"

#include "brg/physics/astro/astro.h"
#include "brg/physics/astro/density_profile/density_profile.h"
#include "brg/physics/astro/density_profile/tNFW_profile.h"
#include "brg/physics/astro/SALTSA/gabdt.h"
#include "brg/physics/phase.hpp"
#include "brg/physics/units/unit_obj.h"

namespace brgastro
{

/** Class Forward Declarations **/
#if (1)
// Forward declare stripping_orbit_segment, whose header has to be included at the bottom of this file
// for the enum used here to work correctly
class stripping_orbit_segment;

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

	/**
	 * An enum of the possible types of interpolation the class will use
	 * to estimate intermediate points on the satellite's orbit. If not
	 * set explicitly, it will default to UNSET, which chooses the
	 * interpolation type from the others. It will usually choose SPLINE
	 * interpolation, but for a large enough number of known points, it
	 * will choose LINEAR interpolation, and for an even larger number it
	 * will choose LOWER.
	 */
	enum allowed_interpolation_type {
		LOWER, //!< No interpolation, use nearest known point at a lower time.
		UPPER, //!< No interpolation, use nearest known point at a higher time.
		LINEAR,//!< Linear interpolation
		SPLINE,//!< Cubic spline interpolation
		UNSET  //!< Automatically determine optimal interpolation type
	}; // end enum allowed_interpolation_type

private:
#if (1)

	// Default integration parameters
#if(1)
	/** Default number of steps for which stripping is calculated. Implemented in
	 *  brg_orbit.cpp. */
	static size_t _default_spline_resolution_;

	/// Default interpolation method. Implemented in brg_orbit.cpp.
	static allowed_interpolation_type _default_interpolation_type_;

	//@{
	/**
	 * Variable step length tweaking: Time step length is proportional to (v_0/v)^(step_length_power)
	 * This gives smaller steps when the satellite is moving faster.
	 * If you want to turn off adaptive step size, set step_length_power to 0
	 * Alternatively, set step_length_power to 1 for even steps in position
	 */
	static BRG_VELOCITY _default_v_0_; /// "Pivot" velocity
	static BRG_DISTANCE _default_r_0_; /// "Pivot" radial distance
	static double _default_step_length_power_; /// How strongly variable step length is implemented
	static double _default_step_factor_max_; /// Maximum allowed value of (v_0/v)^(step_length_power)
	static double _default_step_factor_min_; /// Minimum allowed value of (v_0/v)^(step_length_power)
	//@}
#endif

	// Integration parameters
#if(1)
	/// Estimated number of steps for which stripping is calculated. The actual number of steps
	/// will depend on the adaptive step size calculations, and will usually be larger than this.
	size_t _base_resolution_;

	/// Interpolation method for this orbit
	allowed_interpolation_type _interpolation_type_;

	//@{
	/**
	 * Variable step length tweaking for this orbit. Time step length is proportional to
	 * (v_0/v)^(step_length_power).
     * This gives smaller steps when the satellite is moving faster.
	 * If you want to turn off adaptive step size, set step_length_power to 0
	 * Alternatively, set step_length_power to 1 for even steps in position
	 */
	BRG_VELOCITY _v_0_; /// "Pivot" velocity
	BRG_DISTANCE _r_0_; /// "Pivot" radial distance
	double _step_length_power_; /// How strongly variable step length is implemented
	double _step_factor_max_; /// Maximum allowed value of (v_0/v)^(step_length_power)
	double _step_factor_min_; /// Minimum allowed value of (v_0/v)^(step_length_power)
	//@}
#endif

	// Default tuning parameters
#if(1)
	//@{
	/// Default tuning parameters, for how strong stripping and shocking are and when shocking is active.
	static double _default_tidal_stripping_amplification_; ///< Amplifies tidal stripping by this factor
	static double _default_tidal_stripping_deceleration_; /**< If positive, increase tidal stripping near pericentre,
														        if negative, decrease near pericentre */
	static double _default_tidal_stripping_radialness_; /**< How much tidal stripping depends on full velocity
	                                                         versus tangential velocity. Larger value of this
	                                                         increases stripping of more radial orbits preferentially */
	static double _default_tidal_shocking_amplification_; ///< Amplifies tidal heating by this factor
	static double _default_tidal_shocking_persistance_; ///< How long shocking is active for
	static double _default_tidal_shocking_power_; ///< Affects interplay of stripping and satellite halo profile
	//@}
#endif

	// Tuning parameters
#if(1)
	//@{
	/// Orbit-specific tuning parameters, for how strong stripping and shocking are and when shocking is active.
	double _tidal_stripping_amplification_; /// Amplifies tidal stripping by this factor
	double _tidal_stripping_deceleration_; /** If positive, increase tidal stripping near pericentre,
	  	  	  	  	  	  	  	  	  	       if negative, decrease near pericentre */
	double _tidal_stripping_radialness_; /** How much tidal stripping depends on full velocity
	                                         versus tangential velocity. Larger value of this
	                                         increases stripping of more radial orbits preferentially */
	double _tidal_shocking_amplification_; /// Amplifies tidal heating by this factor
	double _tidal_shocking_persistance_; /// How long shocking is active for
	double _tidal_shocking_power_; /// Affects interplay of stripping and satellite halo profile
	//@}
#endif

	// Global info for the orbit
#if(1)

	mutable size_t _num_segments_; ///< Number of segments the orbit is split up into.

	//@{
	/// The values of t_min/max based only on the points passed to the orbit.
	BRG_TIME _t_min_natural_value_, _t_max_natural_value_;
	//@}

	//@{
	/// The values of t_min/max the user has told us to use.
	BRG_TIME _t_min_override_value_, _t_max_override_value_;
	//@}

	//@{
	/// Whether or not the user has told us to use specific t_min/max values.
	bool _override_t_min_, _override_t_max_;
	//@}

#endif

	// Data for output info for satellite and host_ptr
#if(1)
	//@{
	/// Unit conversion factors for the satellite and host profiles' parameters.
	std::vector< double > _satellite_parameter_unitconvs_,
			_host_parameter_unitconvs_;
	//@}

	//@{
	/// Which of the satellite and host profiles' parameters will be output.
	std::vector< bool > _satellite_output_parameters_,
			_host_output_parameters_;
	//@}
#endif

	// Lists of points on the orbit and related info
#if(1)

	//@{
	/// Points on the orbit which the user has told us about.
	std::vector< std::pair< double, double > > _x_points_,
			_y_points_, _z_points_, _test_mass_points_;
	std::vector< std::pair< double, double > > _vx_points_,
			_vy_points_, _vz_points_;
	std::vector< double > _vx_unknown_points_,
			_vy_unknown_points_, _vz_unknown_points_,
			_t_points_, _host_param_t_points_;
	std::vector< std::pair< double, std::vector< BRG_UNITS > > > _host_parameter_points_;
	//@}

	/// Any discontinuity times the user has told us about.
	std::vector< BRG_TIME > _discontinuity_times_;

	/// Discontinuity times, cleaned of duplicates and out-of-bounds values, then sorted
	mutable std::vector< BRG_TIME > _cleaned_discontinuity_times_;

	/// Interpolator for calculated retained mass fraction, so we can estimate it at any time.
	mutable brgastro::interpolator _m_ret_interpolator_;

	/// Interpolator for calculated retained virial mass fraction, so we can estimate it at any time.
	mutable brgastro::interpolator _m_vir_ret_interpolator_;

	/// Interpolator for loaded comparison retained mass fraction, so we can estimate it at any time.
	brgastro::interpolator _test_mass_interpolator_;

	/// Interpolator for loaded comparison retained mass fraction error, so we can estimate it at any time.
	brgastro::interpolator _test_mass_error_interpolator_;

	/// List of the final retained mass fractions for each orbit segment.
	mutable std::vector< double > _final_frac_m_ret_list_;

	/// List of the final retained virial mass fractions for each orbit segment.
	mutable std::vector< double > _final_frac_m_vir_ret_list_;

	/// The set of orbit segments that are calculated in turn.
	mutable std::vector< brgastro::stripping_orbit_segment > _orbit_segments_;
#endif

	// Host and satellite pointers and info
#if(1)
	/// Pointer to the initial host halo. It's const, so there's no worry about it being modified.
	const density_profile *_init_host_ptr_;
	/// Pointer to the initial satellite halo. It's const, so there's no worry about it being modified.
	const density_profile *_init_satellite_ptr_;

	/// Whether or not to record detailed orbit information for output.
	mutable bool _record_full_data_;

	/// Whether or not multiple host param points have been given to us.
	bool _host_is_evolving_;

	/// Whether or not the user has loaded a host halo pointer.
	bool _host_loaded_;

	/// Whether or not the user has loaded a satellite halo pointer.
	bool _satellite_loaded_;

	/// Whether or not the host pointer points to a member tNFW_profile
	bool _using_private_init_host_;

	/// Whether or not the satellite pointer points to a member tNFW_profile
	bool _using_private_init_satellite_;

	/** A member tNFW_profile for the host, in case the user doesn't want to
	 *  define one externally. */
	tNFW_profile _private_tNFW_init_host_;

	/** A member tNFW_profile for the host, in case the user doesn't want to
	 *  define one externally. */
	tNFW_profile _private_tNFW_init_satellite_;
#endif

	// Calculation output data
#if(1)
	/// Whether or not calculation has completed (successfully or not).
	mutable bool _calculated_;

	/// Whether or not an error of any sort occurred while calculating stripping.
	mutable bool _bad_result_;

	/** Whether or not an error occurred in a circumstance that would likely cause the satellite
	 *  to have been disrupted. */
	mutable bool _likely_disrupted_;

	/// The final segment for which stripping was successfully calculated.
	mutable std::vector< brgastro::stripping_orbit_segment >::iterator _final_good_segment_;
#endif

	// Private methods
#if(1)
	/**
	 * Load an orbit segment with all necessary parameters and initial data for calculation.
	 *
	 * @param segment The segment to be set up
	 * @param segment_init_satellite The satellite's halo profile at the beginning of this segment
	 * @param segment_init_host The host's halo profile at the beginning of this segment
	 * @param resolution The resolution for this segment
	 */
	void _pass_parameters_to_segment(
			brgastro::stripping_orbit_segment & segment,
			brgastro::density_profile *segment_init_satellite=NULL,
			brgastro::density_profile *segment_init_host=NULL,
			size_t resolution=0) const;

	/**
	 * Calculate if necessary, then return an iterator to the final good segment.
	 *
	 * @return Iterator pointing to the final good segment. On failure, it will be _orbit_segments_.end()
	 */
	const std::vector< brgastro::stripping_orbit_segment >::iterator _final_good_segment() const;

#endif // Private methods

#endif // Private members and methods

public:
#if (1)

	// Swap functions
	/**
	 * Explicitly-defined swap method.
	 *
	 * @param other Other stripping_orbit to be swapped with
	 */
	void swap(stripping_orbit &other);
	/**
	 * Explicitly-defined swap function for stripping orbits.
	 * @param same One of the orbits to be swapped.
	 * @param other The other orbit to be swapped.
	 */
	friend void swap(stripping_orbit &same, stripping_orbit &other) {same.swap(other);}

	// Constructors, destructor, and related functions

	/// Default constructor
	stripping_orbit();

	/**
	 * Constructor with initial haloes and optionally a specific resolution. All this needs
	 * is a set of points on the orbit, and then it can calculate.
	 *
	 * @param init_host Initial host profile.
	 * @param init_satellite Initial satellite profile.
	 * @param base_resolution Base resolution for stripping calculate, before adaptive step size is applied.
	 */
	stripping_orbit( const density_profile *init_host, const density_profile *init_satellite,
			const int base_resolution = 200 );
	/**
	 * Copy constructor. This makes sure all pointers and iterators are copied properly.
	 * @param other_orbit_spline Other stripping_orbit to be copied.
	 */
	stripping_orbit( const stripping_orbit &other_orbit_spline ); // Copy constructor

	/**
	 * Assignment operator
	 *
	 * @param other Orbit to be copied.
	 * @return Reference to self.
	 */
	stripping_orbit & operator=( stripping_orbit other );

	/**
	 * Clone function, which creates a new copy of this.
	 * @return Copy of this, allocated with new. Make sure to delete when done with it.
	 */
	stripping_orbit *stripping_orbit_clone();

	/// Standard virtual destructor
	virtual ~stripping_orbit();

	// Setting default integration parameters
#if(1)
	//@{
	/**
	 * Functions to set default integration parameters. The one-argument version only
	 * affects the default, but can be called without an instance of this class existing.
	 *
	 * The multiple-argument version can override the current value as well if
	 * override_current==true.
	 */

	static void set_default_resolution( const size_t new_default_spline_resolution);
	void set_default_resolution( const size_t new_default_spline_resolution,
			const bool override_current,
			const bool silent=false );
	static void set_default_interpolation_type(
			const allowed_interpolation_type new_default_interpolation_type);
	void set_default_interpolation_type(
			const allowed_interpolation_type new_default_interpolation_type,
			const bool override_current,
			const bool silent=false );
	static void set_default_v_0( const double new_default_v_0);
	void set_default_v_0( const double new_default_v_0,
			const bool override_current,
			const bool silent=false );
	static void set_default_r_0( const double new_default_r_0);
	void set_default_r_0( const double new_default_r_0,
			const bool override_current,
			const bool silent=false );
	static void set_default_step_length_power( const double new_default_step_length_power);
	void set_default_step_length_power( const double new_default_step_length_power,
			const bool override_current,
			const bool silent=false );
	static void set_default_step_factor_max( const double new_default_step_factor_max);
	void set_default_step_factor_max( const double new_default_step_factor_max,
			const bool override_current,
			const bool silent=false );
	static void set_default_step_factor_min( const double new_default_step_factor_min);
	void set_default_step_factor_min( const double new_default_step_factor_min,
			const bool override_current,
			const bool silent=false );
	//@}
#endif

	// Setting default tuning parameters
#if(1)
	//@{
	/**
	 * Functions to set default tuning parameters. The one-argument version only
	 * affects the default, but can be called without an instance of this class existing.
	 *
	 * The multiple-argument version can override the current value as well if
	 * override_current==true.
	 */
	static void set_default_tidal_stripping_amplification(
			const double new_default_tidal_stripping_amplification);
	void set_default_tidal_stripping_amplification(
			const double new_default_tidal_stripping_amplification,
			const bool override_current,
			const bool silent=false );
	static void set_default_tidal_stripping_deceleration(
			const double new_default_tidal_stripping_deceleration);
	void set_default_tidal_stripping_deceleration(
			const double new_default_tidal_stripping_deceleration,
			const bool override_current,
			const bool silent=false );
	static void set_default_tidal_stripping_radialness(
			const double new_default_tidal_stripping_radialness);
	void set_default_tidal_stripping_radialness(
			const double new_default_tidal_stripping_radialness,
			const bool override_current,
			const bool silent=false );
	static void set_default_tidal_shocking_amplification(
			const double new_default_tidal_shocking_amplification);
	void set_default_tidal_shocking_amplification(
			const double new_default_tidal_shocking_amplification,
			const bool override_current,
			const bool silent=false );
	static void set_default_tidal_shocking_persistance(
			const double new_default_tidal_shocking_persistance);
	void set_default_tidal_shocking_persistance(
			const double new_default_tidal_shocking_persistance,
			const bool override_current,
			const bool silent=false );
	static void set_default_tidal_shocking_power(
			const double new_default_tidal_shocking_power);
	void set_default_tidal_shocking_power(
			const double new_default_tidal_shocking_power,
			const bool override_current,
			const bool silent=false );
	//@}
#endif

	// Setting integration parameters
#if(1)
	/**
	 * Set the base resolution for this orbit (before adaptive step size adjustments)
	 *
	 * @param new_spline_resolution The new resolution.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_resolution( const size_t new_spline_resolution,
			const bool silent=false );
	/**
	 * Set the type of interpolation to be used for all interpolators, using the
	 * allowed_interpolation_type enum.
	 *
	 * @param new_type The new interpolation type.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_interpolation_type( const allowed_interpolation_type new_type,
			const bool silent=false );
	/**
	 * Set the "pivot" velocity for adaptive step size calculations. If the satellite
	 * is travelling at this velocity and at the pivot radial distance, no step size
	 * correction will be imposed. The step size will be the lower of
	 * (v_0/v)^(step_length_power) and (r/r_0)^(step_length_power), within the bounds
	 * of the min and max.
	 *
	 * @param new_v_0 The new "pivot" velocity.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_v_0( const double new_v_0,
			const bool silent=false );
	/**
	 * Set the "pivot" radial distance for adaptive step size calculations. If the satellite
	 * is travelling at this radial distance and at the pivot velocity, no step size
	 * correction will be imposed. The step size will be the lower of
	 * (v_0/v)^(step_length_power) and (r/r_0)^(step_length_power), within the bounds
	 * of the min and max.
	 *
	 * @param new_r_0 The new "pivot" radial distance.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_r_0( const double new_r_0,
			const bool silent=false );
	/**
	 * Set how powerful the adaptive step size correction is.
	 * The step size will be the lower of (v_0/v)^(step_length_power) and
	 * (r/r_0)^(step_length_power), within the bounds of the min and max.
	 *
	 * @param new_step_length_power The new step length power.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_step_length_power( const double new_step_length_power,
			const bool silent=false );
	/**
	 * Set the maximum allowed step size factor. eg. if this is 10, the
	 * step size will be at most 10 times the size it would be with no
	 * correction.
	 *
	 * @param new_step_factor_max The new maximum step length factor.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_step_factor_max( const double new_step_factor_max,
			const bool silent=false );
	/**
	 * Set the minimum allowed step size factor. eg. if this is 0.1, the
	 * step size will be at minimum 0.1 times the size it would be with no
	 * correction.
	 *
	 * @param new_step_factor_min The new minimum step length factor.
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_step_factor_min( const double new_step_factor_min,
			const bool silent=false );
#endif

	// Setting tuning parameters
#if(1)
	// Tuning parameters, for how strong stripping and shocking are and when shocking is active
	/**
	 * Set the new value for tidal stripping amplification (A_s in the paper).
	 *
	 * @param new_tidal_stripping_amplification
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_stripping_amplification( const double new_tidal_stripping_amplification,
			const bool silent=false );
	/**
	 * Set the new value for tidal stripping deceleration (alpha_s in the paper).
	 *
	 * @param new_tidal_stripping_deceleration
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_stripping_deceleration( const double new_tidal_stripping_deceleration,
			const bool silent=false );
	/**
	 * Set the new value for tidal stripping radialness (alpha_v in the paper).
	 *
	 * @param new_tidal_stripping_deceleration
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_stripping_radialness( const double new_tidal_stripping_radialness,
			const bool silent=false );
	/**
	 * Set the new value for tidal shocking amplification (A_h in the paper).
	 *
	 * @param new_tidal_shocking_amplification
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_shocking_amplification( const double new_tidal_shocking_amplification,
			const bool silent=false );
	/**
	 * Set the new value for tidal shocking persistance (fixed to 1 in the paper).
	 *
	 * @param new_tidal_shocking_persistance
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_shocking_persistance( const double new_tidal_shocking_persistance,
			const bool silent=false );
	/**
	 * Set the new value for tidal shocking power (alpha_h in the paper).
	 *
	 * @param new_tidal_shocking_power
	 * @param silent Whether or not to suppress error messages.
	 */
	void set_tidal_shocking_power( const double new_tidal_shocking_power,
			const bool silent=false );
#endif

	// Resetting integration parameters
#if(1)
	//@{
	/// Reset integration parameters to their default values
	void reset_resolution();
	void reset_interpolation_type();
	void reset_v_0();
	void reset_r_0();
	void reset_step_length_power();
	void reset_step_factor_max();
	void reset_step_factor_min();
	//@}
#endif

	// Resetting tuning parameters
#if(1)
	//@{
	/// Reset tuning parameters to their default values
	void reset_tidal_stripping_amplification();
	void reset_tidal_stripping_deceleration();
	void reset_tidal_stripping_radialness();
	void reset_tidal_shocking_amplification();
	void reset_tidal_shocking_persistance();
	void reset_tidal_shocking_power();
	//@}
#endif

	// Adding data points and clearing those vectors
#if(1)

	/**
	 * Tell the orbit about a point along the satellite's path relative to the host, where
	 * all position and velocity values are known. Optionally also tell it about a comparison
	 * retained mass fraction at this point and the error on that value. This version checks
	 * to ensure that no duplicate times are added. This means that adding a total of N
	 * points will take O(N^2) time (albeit with a low coefficient). If you're sure that there
	 * are no duplicate times, you can use the force_add_point version. An exception will be
	 * thrown if a duplicate time is added here.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input position values which you have in kpc, use
	 * x*unitconv::kpctom.
	 *
	 * @param x x at this point in the orbit.
	 * @param y y at this point in the orbit.
	 * @param z z at this point in the orbit.
	 * @param vx vx at this point in the orbit.
	 * @param vy vy at this point in the orbit.
	 * @param vz vz at this point in the orbit.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param test_mass Optional: A comparison retained mass fraction at this point in the orbit.
	 * @param test_mass_error Optional: Error on the comparison retained mass fraction.
	 */
	void add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
			CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
			const double test_mass = 1, const double test_mass_error = 1 );
	/**
	 * Tell the orbit about a point along the satellite's path relative to the host, where
	 * all position and velocity values are known. Optionally also tell it about a comparison
	 * retained mass fraction at this point and the error on that value. This version does
	 * not check for duplicate times, so only use it if you're sure no times will be duplicates.
	 * If a duplicate time IS added here, it will be caught later, but the exception at that point
	 * will be buried within the interpolator and won't help you catch where the issue is.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input position values which you have in kpc, use
	 * x*unitconv::kpctom.
	 *
	 * @param x x at this point in the orbit.
	 * @param y y at this point in the orbit.
	 * @param z z at this point in the orbit.
	 * @param vx vx at this point in the orbit.
	 * @param vy vy at this point in the orbit.
	 * @param vz vz at this point in the orbit.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param test_mass Optional: A comparison retained mass fraction at this point in the orbit.
	 * @param test_mass_error Optional: Error on the comparison retained mass fraction.
	 */
	void force_add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_VELOCITY_REF vx,
			CONST_BRG_VELOCITY_REF vy, CONST_BRG_VELOCITY_REF vz, CONST_BRG_TIME_REF t,
			const double test_mass = 1, const double test_mass_error = 1 );
	/**
	 * Tell the orbit about a point along the satellite's path relative to the host, where
	 * Position is known but not velocity. The class will attempt to calculate unknown velocity
	 * points through differentiating the position interpolators, but this is time-consuming, so
	 * only use this version when necessary. Optionally also tell it about a comparison
	 * retained mass fraction at this point and the error on that value. This version checks
	 * to ensure that no duplicate times are added. This means that adding a total of N
	 * points will take O(N^2) time (albeit with a low coefficient). If you're sure that there
	 * are no duplicate times, you can use the force_add_point version. An exception will be
	 * thrown if a duplicate time is added here.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input position values which you have in kpc, use
	 * x*unitconv::kpctom.
	 *
	 * @param x x at this point in the orbit.
	 * @param y y at this point in the orbit.
	 * @param z z at this point in the orbit.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param test_mass Optional: A comparison retained mass fraction at this point in the orbit.
	 * @param test_mass_error Optional: Error on the comparison retained mass fraction.
	 */
	void add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
			const double new_test_mass = 1, const double test_mass_error = 1 ); // Only use if v is unknown
	/**
	 * Tell the orbit about a point along the satellite's path relative to the host, where
	 * Position is known but not velocity. The class will attempt to calculate unknown velocity
	 * points through differentiating the position interpolators, but this is time-consuming, so
	 * only use this version when necessary. Optionally also tell it about a comparison
	 * retained mass fraction at this point and the error on that value. This version does
	 * not check for duplicate times, so only use it if you're sure no times will be duplicates.
	 * If a duplicate time IS added here, it will be caught later, but the exception at that point
	 * will be buried within the interpolator and won't help you catch where the issue is.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input position values which you have in kpc, use
	 * x*unitconv::kpctom.
	 *
	 * @param x x at this point in the orbit.
	 * @param y y at this point in the orbit.
	 * @param z z at this point in the orbit.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param test_mass Optional: A comparison retained mass fraction at this point in the orbit.
	 * @param test_mass_error Optional: Error on the comparison retained mass fraction.
	 */
	void force_add_point( CONST_BRG_DISTANCE_REF x, CONST_BRG_DISTANCE_REF y,
			CONST_BRG_DISTANCE_REF z, CONST_BRG_TIME_REF t,
			const double new_test_mass = 1, const double test_mass_error = 1 ); // Only use if v is unknown
	/**
	 * Tell the orbit about the time of a discontinuity. For instance, this could be where the
	 * host group changes.
	 *
	 * @param t The time of a discontinuity in seconds.
	 *          Use the function SALTSA::zft(z) if only redshift is known.
	 */
	void add_discontinuity_time( CONST_BRG_TIME_REF t ); // Splits into segments to be calculated individually
	/**
	 * Tell the orbit about the state of the host halo at a given time. This version checks
	 * to ensure that no duplicate times are added (only among host parameter point times).
	 * This means that adding a total of N points will take O(N^2) time (albeit with a low coefficient).
	 * If you're sure that there are no duplicate times, you can use the force_add_point version.
	 * An exception will be thrown if a duplicate time is added here.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input mass values which you have in 10^10 Msun, use
	 * M*unitconv::ttMsuntom.
	 *
	 * @param parameters A vector of parameters which define the state of the host group at this time.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param silent Whether or not to suppress error messages.
	 */
	void add_host_parameter_point( const std::vector< BRG_UNITS > &parameters, CONST_BRG_TIME_REF t,
			const bool silent = false );
	/**
	 * Tell the orbit about the state of the host halo at a given time. This version does
	 * not check for duplicate times, so only use it if you're sure no times will be duplicates.
	 * If a duplicate time IS added here, it will be caught later, but the exception at that point
	 * will be buried within the interpolator and won't help you catch where the issue is.
	 *
	 * Make sure all values are added in SI units. If you don't have them in these units, use
	 * the unitconv namespace. For instance, to input mass values which you have in 10^10 Msun, use
	 * M*unitconv::ttMsuntom.
	 *
	 * @param parameters A vector of parameters which define the state of the host group at this time.
	 * @param t The time of this point. If only redshift is available, use the function SALTSA::zft(z).
	 * @param silent Whether or not to suppress error messages.
	 */
	void force_add_host_parameter_point( const std::vector< BRG_UNITS > &parameters, CONST_BRG_TIME_REF t,
			const bool silent = false );
	/**
	 * Clears all phase-space points that the orbit has been told about.
	 */
	void clear_points();
	/**
	 * Clears all discontinuity times that the orbit has been told about.
	 */
	void clear_discontinuity_times();
	/**
	 * Clears all host parameter points that the orbit has been told about.
	 */
	void clear_host_parameter_points();

#endif // Adding data points and clearing those vectors

	// Setting and clearing _init satellite/host_ptr
#if(1)
	/**
	 * Set a profile for the initial state of the satellite halo.
	 *
	 * @param new_init_satellite Const pointer to a derived class of density_profile.
	 */
	void set_init_satellite( const density_profile *new_init_satellite );
	/**
	 * Set a profile for the initial state of the host halo.
	 *
	 * @param new_init_host Const pointer to a derived class of density_profile.
	 */
	void set_init_host( const density_profile *new_init_host );
	/**
	 * Tell the orbit to generate its own tNFW halo to be used for the initial state of the satellite.
	 *
	 * @param new_init_mvir0 Mvir0 value for the profile in kg (You can simply use virial mass here).
	 * @param z Redshift of the profile; defaults to 0
	 * @param new_init_c Concentration of the profile; will be determined through a fitting function with
	 *                   mass if left to default.
	 * @param new_init_tau Truncation parameter tau of the profile. Will be 2*c if left to default
	 */
	void set_tNFW_init_satellite( CONST_BRG_MASS_REF new_init_mvir0,
			const double z = 0, const double new_init_c = -1,
			const double new_init_tau = -1 );
	/**
	 * Tell the orbit to generate its own tNFW halo to be used for the initial state of the host.
	 *
	 * @param new_init_mvir0 Mvir0 value for the profile in kg (You can simply use virial mass here).
	 * @param z Redshift of the profile; defaults to 0
	 * @param new_init_c Concentration of the profile; will be determined through a fitting function with
	 *                   mass if left to default.
	 * @param new_init_tau Truncation parameter tau of the profile. Will be 2*c if left to default
	 */
	void set_tNFW_init_host( CONST_BRG_MASS_REF new_mvir0, const double z = 0,
			const double new_c = -1, const double new_tau = -1 );
	/**
	 * Clears initial satellite profile. Not actually necessary to do at any point, as setting to
	 * something else will get the job done, but this is here in case someone wants to use it.
	 */
	void clear_init_satellite();
	/**
	 * Clears initial host profile. Not actually necessary to do at any point, as setting to
	 * something else will get the job done, but this is here in case someone wants to use it.
	 */
	void clear_init_host();

#endif // Setting and clearing init satellite/host_ptr

	// Setting and resetting t_min_/max
#if(1)

	/**
	 * Override the default minimum time for stripping calculation. For instance, use this if
	 * you only want to calculate stripping for a subsection of the orbit.
	 *
	 * @param new_t_min The new minimum time.
	 */
	void set_t_min( CONST_BRG_TIME_REF  new_t_min );
	/**
	 * Override the default maximum time for stripping calculation. For instance, use this if
	 * you only want to calculate stripping for a subsection of the orbit.
	 *
	 * @param new_t_min The new maximum time.
	 */
	void set_t_max( CONST_BRG_TIME_REF  new_t_max );
	/**
	 * Reset the minimum time to the default value (which will be the lowest time point passed
	 * to the orbit).
	 */
	void reset_t_min();
	/**
	 * Reset the maximum time to the default value (which will be the highest time point passed
	 * to the orbit).
	 */
	void reset_t_max();

#endif // Setting and resetting _t_min_/max

	// Functions for determining how calc() will be called
	/**
	 * Tell the orbit whether or not you want it to store full data for the state of the orbit.
	 * If you plan to print out data, set this to true before calculation to save time. The orbit
	 * defaults to not storing this data to save memory.
	 *
	 * @param new_record_full_data Whether or not to save full data on the orbit.
	 */
	void set_record_full_data( const bool new_record_full_data ) const;

	// Global clearing functions
	/**
	 * Clear everything about the orbit, making it as good as new.
	 */
	void clear();
	/**
	 * Clears all calculated results from the orbit.
	 */
	void clear_calcs() const;

	/**
	 * Tell the orbit to calculate stripping now.
	 *
	 * @param silent Whether or not to suppress error messages.
	 */
	void calc( const bool silent = false ) const; // Using conceptual const-ness

	// Output-modifying functions
#if(1)

	/**
	 * Tell the orbit which of the satellite halo's defining parameters you want to be
	 * output.
	 *
	 * @param satellite_output_parameters Vector of bools, saying whether each parameter should be
	 *                                    output or not.
	 */
	void set_satellite_output_parameters(
			const std::vector< bool > &satellite_output_parameters );
	/**
	 * Tell the orbit what unit conversions to use in outputting the satellite halo's
	 * defining parameters.
	 *
	 * @param num_parameters Number of parameters which define the satellite profile.
	 * @param satellite_unitconvs vector of unit conversion factors to use. Use the unitconv namespace
	 *                            for these, with the form unitconv::ttMsuntokg to output mass in
	 *                            10^10 Msun, for instance (units you want first).
	 */
	void set_satellite_parameter_unitconvs(
			const std::vector< double > &satellite_unitconvs );
	/**
	 * Tell the orbit which of the host halo's defining parameters you want to be
	 * output.
	 *
	 * @param host_output_parameters Vector of bools, saying whether each parameter should be
	 *                                    output or not.
	 */
	void set_host_output_parameters( const std::vector< bool > &host_output_parameters );
	/**
	 * Tell the orbit what unit conversions to use in outputting the host halo's
	 * defining parameters.
	 *
	 * @param host_unitconvs vector of unit conversion factors to use. Use the unitconv namespace
	 *                            for these, with the form unitconv::ttMsuntokg to output mass in
	 *                            10^10 Msun, for instance (units you want first).
	 */
	void set_host_parameter_unitconvs( const std::vector< double > &host_unitconvs );
	/**
	 * Clear the vector of which satellite parameters to output. This will result in the default
	 * behaviour of all being output.
	 */
	void clear_satellite_output_parameters();
	/**
	 * Clear the vector of what unit conversions to use for the satellite's parameters. This will
	 * result in all parameters being output in SI units (or as they are if they're unitless).
	 */
	void clear_satellite_parameter_unitconvs();
	/**
	 * Clear the vector of which host parameters to output. This will result in the default
	 * behaviour of all being output.
	 */
	void clear_host_output_parameters();
	/**
	 * Clear the vector of what unit conversions to use for the host parameters. This will
	 * result in all parameters being output in SI units (or as they are if they're unitless).
	 */
	void clear_host_parameter_unitconvs();
#endif

	/**
	 * Print data for the full orbit into the given stream. What data will be output and what
	 * units it uses are determined by the set_output_parameters and set_output_parameter_unitconvs
	 * functions.
	 *
	 * @param out Out stream to print the orbit data into (for instance, an open ofstream).
	 */
	void print_full_data( std::ostream *out ) const;
	/**
	 * Print data for the a specific segment into the given stream. What data will be output and what
	 * units it uses are determined by the set_output_parameters and set_output_parameter_unitconvs
	 * functions.
	 *
	 * @param out Out stream to print the orbit data into (for instance, an open ofstream).
	 * @param segment_number Index for the segment you want to be output
	 */
	void print_segment_data( std::ostream *out,
			const int segment_number ) const;

	// Accessors to private data
#if(1)

	// Default integration parameters
#if(1)
	//@{
	/// Accessors to default integration parameters.
	static const size_t default_spline_resolution() {return _default_spline_resolution_;}
	static const allowed_interpolation_type default_interpolation_type()
		{return _default_interpolation_type_;}
	static BRG_VELOCITY  default_v_0() {return _default_v_0_;}
	static BRG_DISTANCE  default_r_0() {return _default_r_0_;}
	static const double  default_step_length_power() {return _default_step_length_power_;}
	static const double  default_step_factor_max() {return _default_step_factor_max_;}
	static const double  default_step_factor_min() {return _default_step_factor_min_;}
	//@}
#endif

	// Default tuning parameters
#if(1)
	//@{
	/// Accessors to default tuning parameters.
	static const double  default_tidal_stripping_amplification() {return _default_tidal_stripping_amplification_;}
	static const double  default_tidal_stripping_deceleration() {return _default_tidal_stripping_deceleration_;}
	static const double  default_tidal_stripping_radialness() {return _default_tidal_stripping_radialness_;}
	static const double  default_tidal_shocking_amplification() {return _default_tidal_shocking_amplification_;}
	static const double  default_tidal_shocking_persistance() {return _default_tidal_shocking_persistance_;}
	static const double  default_tidal_shocking_power() {return _default_tidal_shocking_power_;}
	//@}
#endif


	// Integration parameters
#if(1)
	//@{
	/// Accessors to integration parameters for this orbit.
	const size_t spline_resolution() const {return _base_resolution_;}
	const allowed_interpolation_type interpolation_type()
		{return _interpolation_type_;}
	BRG_VELOCITY v_0() const {return _v_0_;}
	BRG_DISTANCE r_0() const {return _r_0_;}
	const double step_length_power() const {return _step_length_power_;}
	const double step_factor_max() const {return _step_factor_max_;}
	const double step_factor_min() const {return _step_factor_min_;}
	//@}
#endif

	// Tuning parameters
#if(1)
	//@{
	/// Accessors to tuning parameters for this orbit.
	const double tidal_stripping_amplification() const {return _tidal_stripping_amplification_;}
	const double tidal_stripping_deceleration() const {return _tidal_stripping_deceleration_;}
	const double tidal_stripping_radialness() const {return _tidal_stripping_radialness_;}
	const double tidal_shocking_amplification() const {return _tidal_shocking_amplification_;}
	const double tidal_shocking_persistance() const {return _tidal_shocking_persistance_;}
	const double tidal_shocking_power() const {return _tidal_shocking_power_;}
	//@}
#endif

	/// Accessor to the number of segments this orbit has.
	const size_t num_segments() const {return _num_segments_;}

	/// Accessor to what t_min would be if not overridden.
	BRG_TIME t_min_natural_value() const {return _t_min_natural_value_;}
	/// Accessor to what t_max would be if not overridden.
	BRG_TIME t_max_natural_value() const {return _t_max_natural_value_;}
	/// Accessor to the override value of t_min
	BRG_TIME t_min_override_value() const {return _t_min_override_value_;}
	/// Accessor to the override value of t_max
	BRG_TIME t_max_override_value() const {return _t_max_override_value_;}
	/// Accessor to the current value of t_min
	BRG_TIME t_min() const
	{
		if(_override_t_min_)
			return _t_min_override_value_;
		else
			return _t_min_natural_value_;
	}
	/// Accessor to the current value of t_max
	BRG_TIME t_max() const
	{
		if(_override_t_max_)
			return _t_max_override_value_;
		else
			return _t_max_natural_value_;
	}
	/// Accessor to whether or not to use a user-defined value of t_min
	const bool override_t_min() const {return _override_t_min_;}
	/// Accessor to whether or not to use a user-defined value of t_max
	const bool override_t_max() const {return _override_t_max_;}

	/// Accessor to vector of satellite halo parameter unit conversion factors
	std::vector< double > satellite_parameter_unitconvs() const {return _satellite_parameter_unitconvs_;}
	/// Accessor to vector of host halo parameter unit conversion factors
	std::vector< double > host_parameter_unitconvs() const {return _host_parameter_unitconvs_;}
	/// Accessor to vector of which satellite halo parameters to output
	std::vector< bool > satellite_output_parameters() const {return _satellite_output_parameters_;}
	/// Accessor to vector of which host halo parameters to output
	std::vector< bool > host_output_parameters() const {return _host_output_parameters_;}

	//@{
	/// Accessors to stored data on points passed to the orbit.
	std::vector< std::pair< double, double > > x_spline_points() const {return _x_points_;}
	std::vector< std::pair< double, double > > y_spline_points() const {return _y_points_;}
	std::vector< std::pair< double, double > > z_spline_points() const {return _z_points_;}
	std::vector< std::pair< double, double > > vx_spline_points() const {return _vx_points_;}
	std::vector< std::pair< double, double > > vy_spline_points() const {return _vy_points_;}
	std::vector< std::pair< double, double > > vz_spline_points() const {return _vz_points_;}
	std::vector< double > vx_spline_unknown_points() const {return _vx_unknown_points_;}
	std::vector< double > vy_spline_unknown_points() const {return _vy_unknown_points_;}
	std::vector< double > vz_spline_unknown_points() const {return _vz_unknown_points_;}
	std::vector< std::pair< double, double > > test_mass_spline_points() const
			{return _test_mass_points_;}
	std::vector< std::pair< double, std::vector< BRG_UNITS > > > host_parameter_spline_points() const
			{return _host_parameter_points_;}
	//@}

	/// Accessor to the list of discontinuity times
	std::vector< BRG_TIME > discontinuity_times() const {return _discontinuity_times_;}

	/// Accessor to initial satellite halo pointer
	const density_profile * init_satellite_ptr() const {return _init_satellite_ptr_;}
	/// Accessor to initial host halo pointer
	const density_profile * init_host_ptr() const {return _init_host_ptr_;}

	/// Accessor to whether or not to record full data
	const bool record_full_data() const {return _record_full_data_;}
	/// Accessor to whether or not the host is changing over time
	const bool host_is_evolving() const {return _host_is_evolving_;}
	/// Accessor to whether or not the satellite is loaded
	const bool satellite_loaded() const {return _satellite_loaded_;}
	/// Accessor to whether or not the host is loaded
	const bool host_loaded() const {return _host_loaded_;}
	/// Accessor to whether or not the orbit has been calculated
	const bool calculated() const {return _calculated_;}
	/// Accessor to whether or not there was some error in calculation
	const bool bad_result() const {return _bad_result_;}
	/// Accessor to whether or not the orbit is using an internal initial host halo
	const bool using_private_init_satellite() const {return _using_private_init_satellite_;}
	/// Accessor to whether or not the orbit is using an internal initial satellite halo
	const bool using_private_init_host() const {return _using_private_init_host_;}

	/// Accessor to the orbit's private initial satellite halo
	tNFW_profile private_tNFW_init_satellite() const {return _private_tNFW_init_satellite_;}
	/// Accessor to the orbit's private initial host halo
	tNFW_profile private_tNFW_init_host() const {return _private_tNFW_init_host_;}

	/// Accessor to the vector of segments this orbit has
	std::vector<stripping_orbit_segment> orbit_segments() const {return _orbit_segments_;}
#endif

	// Getting resultant data
#if (1)

	// Get final data (returns 1 on failure)
#if(1)
	/**
	 * Get final retained mass, calculating if necessary.
	 *
	 * @param mret Variable to be loaded with final retained mass.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_m_ret( BRG_MASS & m_ret ) const;
	/**
	 * Get final retained virial mass, calculating if necessary.
	 *
	 * @param fmret Variable to be loaded with final retained mass fraction.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_m_vir_ret( BRG_MASS & m_vir_ret ) const;
	/**
	 * Get final retained mass fraction, calculating if necessary.
	 *
	 * @param mret Variable to be loaded with final retained mass.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_frac_m_ret( double & final_frac_m_ret ) const;
	/**
	 * Get final retained virial mass fraction, calculating if necessary.
	 *
	 * @param mret Variable to be loaded with final retained mass.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_frac_m_vir_ret( double & final_frac_m_ret ) const;
	/**
	 * Get the final sum_deltarho value, in a long double variable.
	 * @param final_sum_deltarho Variable to be loaded with final sum_delta_rho.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_sum_deltarho( long double & final_sum_deltarho ) const;
	/**
	 * Get the final sum_deltarho value, in a double variable.
	 * @param final_sum_deltarho Variable to be loaded with final sum_delta_rho.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_sum_deltarho( BRG_UNITS & final_sum_deltarho ) const;
	/**
	 * Get the final sum_gabdt.
	 * @param final_sum_gabdt Variable to be loaded with final sum_gabdt.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_final_sum_gabdt( gabdt & final_sum_gabdt ) const;
	/**
	 * Get the last discontinuity time of the orbit (assumed to be time of last infall)
	 * @param t Variable to be loaded with final infall time.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_last_infall_time( BRG_TIME & t ) const;
	/**
	 * Get a clone of the final satellite profile.
	 * @param final_satellite_clone Pointer to be created with new, pointing to clone (be sure to delete!)
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int clone_final_satellite(
			density_profile * & final_satellite_clone ) const;
	/**
	 * Get a clone of the final host profile.
	 * @param final_host_clone Pointer to be created with new, pointing to clone (be sure to delete!)
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int clone_final_host( density_profile * & final_host_clone ) const;
#endif // Get final data (returns 1 on failure)

	// Get final data (throws exception on failure)
#if(1)
	/**
	 * Get the final retained mass, throwing an exception on failure.
	 *
	 * @return Final retained mass.
	 */
	const BRG_MASS final_m_ret() const;
	/**
	 * Get the final retained mass fraction, throwing an exception on failure.
	 *
	 * @return Final retained mass fraction.
	 */
	const double final_frac_m_ret() const;
	/**
	 * Get the final retained virial mass, throwing an exception on failure.
	 *
	 * @return Final retained mass.
	 */
	const BRG_MASS final_m_vir_ret() const;
	/**
	 * Get the final retained virial mass fraction, throwing an exception on failure.
	 *
	 * @return Final retained mass fraction.
	 */
	const double final_frac_m_vir_ret() const;
	/**
	 * Get the final sum_deltarho, throwing exception on failure.
	 * @return The final sum_deltarho.
	 */
	const BRG_UNITS final_sum_deltarho() const;
	/**
	 * Get the final sum_gabdt, throwing exception on failure.
	 * @return The final sum_gabdt.
	 */
	const gabdt final_sum_gabdt() const;
	/**
	 * Get the last discontinuity time (assumed to be the last time of infall), throwing an
	 * exception on failure.
	 * @return The last discontinuity time.
	 */
	const BRG_TIME last_infall_time() const;
	/**
	 * Clone the final satellite, throwing exception on failure.
	 * @return A clone of the final satellite. Make sure to store this value and delete it later!
	 */
	const density_profile * final_satellite() const;
	/**
	 * Clone the final host, throwing exception on failure.
	 * @return A clone of the final host. Make sure to store this value and delete it later!
	 */
	const density_profile * final_host() const;
#endif // Get final data (throws exception on failure)

	/**
	 * Get whether or not there was an error which likely disrupted the satellite halo.
	 * @return Whether or not the calculation met an error in a circumstance that likely would have
	 *         disrupted the satellite.
	 */
	const bool likely_disrupted() const;

	// Get data at arbitrary time (returns 1 on failure)
	/**
	 * Get an estimate of the retained mass at an arbitrary time. _record_full_data_ must be set to
	 * True to get this; the orbit will be calculated again with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @param mret Variable to be loaded with the retained mass estimate.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_m_ret_at_t( CONST_BRG_TIME_REF  t, BRG_MASS & m_ret) const;
	/**
	 * Get an estimate of the retained mass fraction at an arbitrary time. _record_full_data_ must be set to
	 * True to get this; the orbit will be calculated again with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass fraction.
	 * @param fmret Variable to be loaded with the retained mass fraction estimate.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_frac_m_ret_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret ) const;
	/**
	 * Get an estimate of the retained virial mass at an arbitrary time. _record_full_data_ must be set to
	 * True to get this; the orbit will be calculated again with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @param mret Variable to be loaded with the retained mass estimate.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_m_vir_ret_at_t( CONST_BRG_TIME_REF  t, BRG_MASS & m_ret) const;
	/**
	 * Get an estimate of the retained virial mass fraction at an arbitrary time. _record_full_data_ must be set to
	 * True to get this; the orbit will be calculated again with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @param mret Variable to be loaded with the retained mass estimate.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_frac_m_vir_ret_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret ) const;
	/**
	 * Get an estimate of the retained mass fraction of the comparison orbit at an arbitrary time.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @param fmret Variable to be loaded with the retained mass estimate.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_comp_frac_m_ret_at_t( CONST_BRG_TIME_REF  t, double & m_ret) const;
	/**
	 * Get an estimate of the retained mass fraction error of the comparison orbit at an arbitrary time.
	 *
	 * @param t Time for which to get an estimate of the retained mass fraction error.
	 * @param fmret_err Variable to be loaded with the retained mass estimate fraction error.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_comp_frac_m_ret_error_at_t( CONST_BRG_TIME_REF  t, double & frac_m_ret_error) const;

	// Get data at arbitrary time (throws exception on failure)
	/**
	 * Get the retained mass at an arbitrary time, throwing an exception on failure.
	 * _record_full_data_ must be set to True to get this; the orbit will be calculated again
	 * with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @return Retained mass estimate.
	 */
	const BRG_MASS m_ret_at_t(CONST_BRG_TIME_REF  t) const;
	/**
	 * Get the retained mass fraction at an arbitrary time, throwing an exception on failure.
	 * _record_full_data_ must be set to True to get this; the orbit will be calculated again
	 * with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @return Retained mass fraction estimate.
	 */
	const double frac_m_ret_at_t(CONST_BRG_TIME_REF  t) const;
	/**
	 * Get the retained virial mass at an arbitrary time, throwing an exception on failure.
	 * _record_full_data_ must be set to True to get this; the orbit will be calculated again
	 * with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @return Retained mass estimate.
	 */
	const BRG_MASS m_vir_ret_at_t(CONST_BRG_TIME_REF  t) const;
	/**
	 * Get the retained virial mass fraction at an arbitrary time, throwing an exception on failure.
	 * _record_full_data_ must be set to True to get this; the orbit will be calculated again
	 * with it set if necessary.
	 *
	 * @param t Time for which to get an estimate of the retained mass.
	 * @return Retained mass fraction estimate.
	 */
	const double frac_m_vir_ret_at_t(CONST_BRG_TIME_REF  t) const;
	/**
	 * Get the retained mass fraction of the comparison orbit at an arbitrary time, throwing an exception
	 * on failure.
	 *
	 * @param t Time for which to get an estimate of the retained mass fraction.
	 * @return Retained mass fraction estimate.
	 */
	const double comp_frac_m_ret_at_t(CONST_BRG_TIME_REF  t) const;
	/**
	 * Get the retained mass fraction error of the comparison orbit at an arbitrary time, throwing an
	 * exception on failure.
	 *
	 * @param t Time for which to get an estimate of the retained mass fraction error.
	 * @return Retained mass fraction error estimate.
	 */
	const double comp_frac_m_ret_error_at_t(CONST_BRG_TIME_REF  t) const;

	// Quality of fit.
	/**
	 * Get a quality-of-fit estimate, comparing the calculated fmret to the comparison values and errors.
	 * The value is taken at a number of sample points and normalized. _record_full_data_ must be set to
	 * True to get this; if it isn't the orbit will be recalculated.
	 *
	 * @param Q Variable to store the quality-of-fit estimate.
     * @param use_virial Set to true to compare retained virial mass fraction instead
	 * @param samples Number of sample points along the orbit to compare fmret values at.
	 * @return Int flag - zero for success, otherwise for error.
	 */
	const int get_quality_of_fit( double & Q, const bool use_virial=false, const unsigned int samples=100 );
	/**
	 * Get a quality-of-fit estimate, comparing the calculated fmret to the comparison values and errors.
	 * The value is taken at a number of sample points and normalized. An exception will be thrown on
	 * error here.  _record_full_data_ must be set to True to get this; if it isn't the orbit will
	 * be recalculated.
	 *
     * @param use_virial Set to true to compare retained virial mass fraction instead
	 * @param samples Number of sample points along the orbit to compare fmret values at.
	 * @return The quality-of-fit estimate.
	 */
	const double quality_of_fit( const bool use_virial=false, const unsigned int samples=100 );

#endif // Getting resultant data

#endif

};
// class stripping_orbit

#endif // end class definitions

} // end namespace brgastro

#include "stripping_orbit_segment.h"

#endif //__BRG_ORBIT_H_INCLUDED__

