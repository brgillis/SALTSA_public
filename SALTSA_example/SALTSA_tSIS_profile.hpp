/**********************************************************************\
  @file SALTSA_tSIS_profile.hpp

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

#ifndef _SALTSA_TSIS_PROFILE_HPP_INCLUDED_
#define _SALTSA_TSIS_PROFILE_HPP_INCLUDED_

#include <cstdlib>
#include <string>
#include <exception>

#include "brg/global.h"

#include "brg/physics/astro.h"
#include "brg/physics/density_profile/density_profile.h"
#include "brg/physics/units/unit_obj.h"
#include "brg/math/misc_math.hpp"

namespace SALTSA {

/**
 * Truncated singular isothermal sphere (tSIS) density profile. This is an example of
 * how you can implement a user-defined density profile to be used with SALTSA.
 *
 * We'll define this profile as having a constant velocity dispersion within the viral
 * radius and zero mass density outside the virial radius.
 */
class tSIS_profile: public brgastro::density_profile {
private:
	/**
	 *  This profile can be defined entirely by redshift and velocity dispersion.
	 *  Redshift is inherited, so we just need to define a variable for velocity
	 *  dispersion, which we'll call sigma_v.
	 */
	double _sigma_v_;
public:

	// Constructors and destructors
#if(1)

	/** A simple constructor, allowing default creation of an empty profile or
	 * through defining the initial velocity dispersion and redshift.
	 *
	 * @param init_sigma_v Initial value for velocity dispersion in m/s.
	 * @param init_z Initial value for redshift.
	 */
	tSIS_profile( double init_sigma_v=0, double init_z=0 )
		: redshift_obj(init_z),
		  _sigma_v_(init_sigma_v)
	{
	}

	// Implicit copy constructor is fine here, since we aren't managing a resource.

	/**
	 *  Virtual destructor. We're not managing a resource with this class, so we
	 *  don't need to do anything here. But we make it virtual out of good practice,
	 *  in case someone wants to inherit from this class.
	 */
	virtual ~tSIS_profile()
	{
	}

#endif // Constructors and destructors

	// New functions for this profile
#if(1)

	/**
	 * Function to set the value of sigma_v. Will throw an exception on a negative or bad
	 * value.
	 *
	 * @param new_sigma_v New value for sigma_v.
	 */
	const int set_sigma_v(const double new_sigma_v)
	{
		if((new_sigma_v<0) || (brgastro::isbad(new_sigma_v)))
			throw std::runtime_error("Cannot set sigma_v to a negative value.");

		_sigma_v_ = new_sigma_v;

		// Since we changed here, tell the profile that cached half-mass radii are no longer
		// valid.
		hmvir_cached = false;
		hmtot_cached = false;

		return 0;
	}

	/**
	 * Accessor to _sigma_v_
	 * @return _sigma_v_
	 */
	double sigma_v()
	{
		return _sigma_v_;
	}


#endif

	// Critical overloaded functions
#if(1)
	/**
	 * Function to get the virial mass. If we assume that the velocity dispersion is
	 * equal to the circular orbital velocity at the virial radius (where the
	 * enclosed density is 200 times the critical density), then we get the
	 * form below.
	 *
	 * @return The virial mass in kg
	 */
	BRG_MASS mvir() const // Virial mass
	{
		return std::pow(_sigma_v_,3.) / (10 * brgastro::Gc * H());
	}

	/**
	 * Function to set all critical parameters of the profile through a vector.
	 *
	 * @param new_parameters The parameter vector
	 * @param silent Set to true to silence any errors (not relevent in this implementation)
	 */
	void set_parameters(const std::vector< BRG_UNITS > & new_parameters,
			const bool silent = false )
	{
		if(new_parameters.size() != num_parameters())
			throw std::logic_error("Invalid number of parameters for tSIS profile: Exactly 2 required.");
		set_sigma_v(new_parameters[0]);
		set_z(new_parameters[1]);
	}

	/**
	 * Function to get the density at a given radius. The form of an SIS's density
	 * profile is pretty standard; we just expand and add the truncation here.
	 *
	 * @param r The radiuse we want to get the density at.
	 * @return The density at this radius, in kg/m^3
	 */
	BRG_UNITS dens( CONST_BRG_DISTANCE_REF r ) const // Local density at radius r
	{
		if(r>rvir())
			return 0;
		else
			return std::pow(_sigma_v_,2.)/(2*pi*brgastro::Gc*pow(r,2.));
	}

	/**
	 * Clone function. This will be used when a deep copy needs to be made of a
	 * pointer to a density profile. We use polymorphism here to return a pointer
	 * to a new copy of this. Since this profile inherits from density profile,
	 * anything expected a density profile clone can take the return from this
	 * function.
	 *
	 * @return A density_profile * which points to a new copy of this (Make sure to delete!)
	 */
	density_profile *density_profile_clone() const // Creates a clone of this
	{
		return new tSIS_profile( *this );
	}
#endif // Critical overloaded functions

	/**
	 * All we actually need to define are the four critical functions here, which are
	 * declared as pure virtual in the density_profile class. If we stop here though,
	 * the class will have to integrate certain other needed values and will be very
	 * slow, so we'll implement some analytic forms that we know. We'll also be missing
	 * out on the ability to use some useful functions, such as setting the mass directly
	 * instead of setting the velocity dispersion
	 */

	// Virtual functions we may be using, and so need to overload
#if(1)

	/**
	 * Function to set mvir, in case we want to set it directly instead of
	 * changing sigma_v. Will throw an exception on negative or bad mvir.
	 *
	 * @param new_mvir New mvir value
	 * @param silent Whether or not to silence errors.
	 */
	void set_mvir( CONST_BRG_MASS_REF new_mvir,
			const bool silent=false )
	{
		if((new_mvir<0) || (brgastro::isbad(new_mvir)))
			throw std::runtime_error("Cannot set mvir to a negative value.");

		_sigma_v_ = std::pow(10 * brgastro::Gc * H() * new_mvir, 1./3.);

		// Since we changed here, tell the profile that cached half-mass radii are no longer
		// valid.

		hmvir_cached = false;
		hmtot_cached = false;
	}

	/**
	 * Function to set vvir, since it's easy enough, being identical to sigma_v.
	 *
	 * @param new_vvir
	 * @param silent
	 */
	void set_vvir( CONST_BRG_VELOCITY_REF new_vvir,
			const bool silent=false )
	{
		set_sigma_v(new_vvir);
	}

	/**
	 * Function to get mtot. Here's it's simply equal to mvir.
	 *
	 * @return Total mass of the halo.
	 */
	BRG_MASS mtot() const // Total mass
	{
		return mvir();
	}

	/**
	 * Function to get tidal radius. Here again, it's simply the virial radius.
	 *
	 * @param silent
	 * @return The truncation radius
	 */
	BRG_DISTANCE rt( const bool silent = false ) const // Tidal/truncation radius
	{
		return rvir();
	}

	/**
	 * Function to get the number of parameters which define this profile. Here we
	 * just have 2: Sigma_v and redshift.
	 *
	 * @return The number of parameters
	 */
	size_t num_parameters() const
	{
		return 2;
	}

	/**
	 * Function to get the parameters which define this profile, which are
	 * Sigma_v and redshift.
	 *
	 * @return
	 */
	std::vector< BRG_UNITS > get_parameters( const bool silent = false ) const // Returns a set of parameters for this halo (all the variables needed to define it - must be defined for each child)
	{
		std::vector< BRG_UNITS > parameters;
		parameters.push_back(_sigma_v_);
		parameters.push_back(z());
		return parameters;
	}

	/**
	 * Function to get the names of the parameters which define this profile, which are
	 * Sigma_v and redshift.
	 *
	 * @param silent
	 * @return A vector of strings, each being the name of one of the parameters
	 */
	std::vector< std::string > get_parameter_names( const bool silent = false ) const // Returns a set of names of this halo's parameters
	{
		std::vector< std::string > parameter_names;
		parameter_names.push_back("sigma_v");
		parameter_names.push_back("z");
		return parameter_names;
	}

	/**
	 * Function to truncate the halo to a specific fraction of its mass. Here we'll implement
	 * this through scaling down the virial mass (which happens to equal total mass here),
	 * using the set_mvir() function we defined above.
	 *
	 * @param fraction to truncate to
	 * @param silent
	 */
	void truncate_to_fraction( const double fraction,
			const bool silent = false ) // Adjusts parameters of this class to decrease the mass to fraction of its previous mass. Must be defined for each child
	{
		set_mvir(fraction*mvir());
	}


#endif

	// Other functions it's useful to override for speed
#if(1)

	/**
	 * Function to get the enclosed mass within a given radius. If we don't implement
	 * this here, it will have to be integrated, which takes up a lot of time. Since
	 * for this profile it's quite simple, we implement it here.
	 *
	 * @param r
	 * @param silent
	 * @return
	 */
	BRG_MASS enc_mass( const double r, const bool silent =
			true ) const // Mass enclosed with sphere of radius r
	{
		if(r>rvir())
			return mvir();
		else
			return 2*r*std::pow(_sigma_v_,2.)/brgastro::Gc;
	}

	/**
	 * Function to get the radius of the sphere which contains half the virial mass.
	 * This too will be integrated if not defined, but it's simple in this case.
	 *
	 * @param silent
	 * @return
	 */
	BRG_DISTANCE rhmvir( const bool silent = false ) const // Half-virial-mass radius
	{
		return rvir()/2;
	}

	/**
	 * Function to get the radius of the sphere which contains half the total mass.
	 * This too will be integrated if not defined, but it's simple in this case.
	 *
	 * @param silent
	 * @return
	 */
	BRG_DISTANCE rhmtot( const bool silent = false ) const // Half-total-mass radius
	{
		return rvir()/2;
	}

	/**
	 * Function to get the virial velocity. It's not a complicated calculation to get
	 * it if it isn't defined, but we can still speed it up.
	 *
	 * @return
	 */
	BRG_VELOCITY vvir() const // Orbital velocity at rvir
	{
		return _sigma_v_;
	}

	/**
	 * Function to get the half-virial-mass orbital velocity. It's not a complicated
	 * calculation to get it if it isn't defined, but we can still speed it up.
	 *
	 * @return
	 */
	BRG_VELOCITY vhmvir( const bool silent = false ) const // Orbital velocity at rhmvir
	{
		return _sigma_v_;
	}

	/**
	 * Function to get the half-total-mass orbital velocity. It's not a complicated
	 * calculation to get it if it isn't defined, but we can still speed it up.
	 *
	 * @return
	 */
	BRG_VELOCITY vhmtot( const bool silent = false ) const // Orbital velocity at rhmtot
	{
		return _sigma_v_;
	}

#endif // Other functions it's useful to override for speed
};

} // end namespace SALTSA

#endif // _SALTSA_TSIS_PROFILE_HPP_INCLUDED_
