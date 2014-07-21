/**********************************************************************\
 @file SALTSA_astro.h
 --------------------

 If this header is used, the source file SALTSA_astro.cpp must be included
 and compiled with the project.

 This file includes various classes and functions for astrophysical
 objects and calculations. The file is split into five primary sections:

 -Constant Definitions
 -Class Forward Declarations
 -Static Class Definitions
 -Function Declarations
 -Class Definitions

 These sections are explained in further detail in their respective
 documentation blocks.

 With the exception of the Constant Definitions section, everything
 in this file is declared in the namespace SALTSA.

 \**********************************************************************/

#ifndef __SALTSA_ASTRO_H_INCLUDED__
#define __SALTSA_ASTRO_H_INCLUDED__

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

#include "SALTSA_global.h"

#include "SALTSA_unitconvs.h"
#include "SALTSA_cache.hpp"
#include "SALTSA_functor.hpp"

/** Constant Definitions **/
#if (1)
/**********************************************************************
 This section stores the definitions of various physical and
 astrophysical constants, stored as variables. Key points:

 -H_0 is set to 70 km/s/Mpc. If this is not changed, all results should
 be assumed to be in h_70 units.
 -Cosmological values are based on best-fit WMAP + priors values
 from Hinshaw et al. 2012.

 \**********************************************************************/

// Physical Constants
#if (1)
#ifndef __PI_DEFINED__
#define __PI_DEFINED__
const double pi = 3.14159265358979323846;
#endif

#ifndef __PHYS_VARS_DEFINED__
#define __PHYS_VARS_DEFINED__

const float Gc = 6.67384e-11; // In m^3 s^-2 kg^-1
const double c = unitconv::ctomps;

#endif

#endif // end Physical Constants

// Cosmological Parameters
#if (1)
#ifndef __COSMO_VARS_DEFINED__
#define __COSMO_VARS_DEFINED__

const double H_0 = 70 * unitconv::kmtom / unitconv::stos / unitconv::Mpctom; // So all results will implicitly be in h_70 units

const float Omega_m = 0.288; // WMAP9 + priors
const float Omega_r = 0.000086; // WMAP9 + priors
const float Omega_k = 0; // Assuming flat space
const float Omega_l = 1 - Omega_k - Omega_m - Omega_r;
const float Omega_b = 0.0472; // WMAP9 + priors
const float sigma_8 = 0.830; // WMAP9 + priors
const float n_s = 0.971;     // WMAP9 + priors

#endif

#endif // end Cosmological Parameters

// Default parameters
#if (1)
#ifndef __BRG_DEFAULT_ASTRO_VARS_DEFINED__
#define __BRG_DEFAULT_ASTRO_VARS_DEFINED__

const double default_c = 6; // To help prevent crashes. Warning will be output
const double default_tau_factor = 2;

#endif
#endif

#endif // End Constant Definitions

namespace SALTSA
{

/** Class Forward Declarations **/
#if (1)
/**********************************************************************
 This section forward declares all class defined in this header, allowing
 them to point to each other if necessary.

 \**********************************************************************/

class redshift_obj;
// Concrete base

class density_profile;
// Abstract base
class tNFW_profile;
// Concrete
class point_mass_profile;
// Concrete

class spherical_density_function;
// Concrete
class accel_function;
// Concrete
class solve_rhm_function;
// Concrete

#endif // end class forward declarations

/** Static Class Definitions **/
#if (1)
/**********************************************************************
 This section defines all "static" classes. That is, classes with
 entirely static variables. These classes are used as caches within the
 program to aid functions which use look-up tables to save time. The use
 of static member variables in this manner ensures that the look-up
 tables only need to be loaded once each.

 \**********************************************************************/

class tfa_cache : public brg_cache<tfa_cache>
{
	// "Time from a (scale factor)" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache;

protected:

	const std::string _name_base() const throw()
	{
		return "tfa";
	}

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class tfa_cache

#endif // End static class definitions

/** Function Declarations **/
#if (1)
/**********************************************************************
 This section defines various functions which I find useful for
 astrophysical calculations. All are declared in the namespace SALTSA.

 \**********************************************************************/

// Functions to work between redshift, scale factor, and time (in s, with zero = present day)
const double zfa( const double a );
const double afz( const double z );

const double tfz( const double z );
const double tfa( const double z );
const double zft( const double t );
const double aft( const double t );

// Functions to integrate out distances
const double integrate_ltd( const double z1, const double z2 );
const double integrate_ltd( const double z );
const double integrate_distance( const double z1, const double z2,
		const int mode, const int resolution = 10000 );

// Functions relating to tNFW profiles
inline const double cfm( const double mass, const double z = 0 ) // Concentration from mass relationship, from Neto
{
	return 4.67 * std::pow( mass * unitconv::kgtottMsun / ( 1e4 ), -0.11 );
}
inline const double mftau( const double tau, const double conc ) // Mtot/Mvir from tau
{
	if(tau<=0) return 0;
	double M0oM200 = 1 / ( log( 1 + conc ) - conc / ( 1 + conc ) );
	double result =
			M0oM200 * tau * tau / std::pow( tau * tau + 1, 2 )
					* ( ( tau * tau - 1 ) * log( tau ) + tau * pi
							- ( tau * tau + 1 ) );
	return result;
}
const double taufm( const double mratio, double c, double tau_init = 0, double precision = 0.00001,
		const bool silent = false ); // tau from Mtot/Mvir
inline const double delta_c( const double conc ) // Simple function of concentration used as a step in calculating NFW densities
{
	return ( 200. / 3. ) * std::pow( conc, 3 )
			/ ( log( 1 + conc ) - conc / ( 1 + conc ) );
}

// Function to estimate orbital period from current position and velocity in a density profile
// Note that this is merely an estimate from analogy to calculations in a Keplerian potential
const double period( const density_profile *host, const double &r,
		const double &vr, const double &vt = 0 );

#endif // end function declarations

/** Class Definitions **/
#if (1)
class redshift_obj
{
	/**********************************
	 redshift_obj class
	 ------------------

	 A base class for anything with a redshift

	 Derived classes:

	 sky_obj
	 galaxy
	 group
	 density_profile
	 tNFW_profile
	 point_mass_profile

	 **********************************/
private:
	double _z_, _z_err_;
public:
	// Constructor
	redshift_obj( const double init_z = 0, const double init_z_err = 0 )
	{
		_z_ = init_z;
		_z_err_ = init_z_err;
	}

	// Copy constructor
	// (implicit is fine for us)

	// Virtual destructor
	virtual ~redshift_obj()
	{
	}

	// Set functions
#if (1)
	virtual const int set_z( const double new_z ) // Sets z
	{
		_z_ = new_z;
		return 0;
	}
	virtual const int set_z_err( const double new_z_err ) // Sets z error
	{
		_z_err_ = new_z_err;
		return 0;
	}
#endif

	// Get functions
#if (1)
	virtual const double z() const
	{
		return _z_;
	}
	virtual const double z_err() const
	{
		return _z_err_;
	}
#endif

	// H(z) at either the redshift of the object or a specified redshift, given in units of m/s^2
	virtual const double H( double ztest = -1 ) const;

	// Clone function - Not needed in current implementation
	// virtual redshift_obj *redshift_obj_clone()=0;
};

class density_profile: public virtual redshift_obj
{
	/**********************************
	 density_profile class
	 -------------

	 An abstract class for anything which has
	 a mass density profile. Meant to be
	 overriden by child classes.

	 The methods for this class are sorted
	 into four groups:
	 -Pure virtual functions (must be implemented
	 by any derived classes)
	 -Virtual functions which must be overwritten if they're going
	 to be used (not pure virtual since they may not be needed)
	 -Virtual functions which should be overwritten if possible and
	 if they're going to be used (those that require integration,
	 but may have an analytic solution for certain profiles)
	 -Virtual functions that don't need to be overwritten in most
	 cases and non-virtual functions non-virtual functions

	 Derived classes:
	 tNFW_profile
	 point_mass profile

	 **********************************/

private:
#if (1)
	int _hm_type_;

	mutable double _rhmvir_cache_, _rhmtot_cache_;

	double _v_from_r( double r ) const
	{
		double a;

		a = accel( fabs( r ) );
		if ( a >= 0 )
			return 0;
		return sqrt( -a * fabs( r ) );
	}
#endif

protected:
	mutable bool hmvir_cached, hmtot_cached;

public:
#if (1)

#if (1) // Constructor
	density_profile()
	{
		_hm_type_ = 1;
		_rhmvir_cache_ = 0;
		hmvir_cached = false;
		_rhmtot_cache_ = 0;
		hmtot_cached = false;
		set_z( 0 );
	}
#endif

#if (1) // Destructor
	virtual ~density_profile()
	{
	}
#endif

#if (1) // Pure virtual functions (must be implemented by any derived classes)

	virtual const double mvir() const =0; // Virial mass (exact definition can be chosen per profile)
	virtual const double dens( const double &r ) const =0; // Local density at radius r

	virtual density_profile *density_profile_clone() const =0; // Creates a clone of this
#endif

#if (1) // Virtual functions which must be overwritten if they're going to be used

#if (1) // Set functions - will return 1 if profile doesn't support this method of setting
	// All take values in default unit set (m, s, kg, K, rad, C)
	virtual const int set_mvir( const double &new_mvir, bool silent = false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_mvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_vvir( const double &new_vvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_vvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rvir( const double &new_rvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rs( const double &new_rs, bool silent = false ) // Scale radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rs(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rt( const double &new_rt, bool silent = false ) // Tidal/truncation radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rt(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< double > &parameters, bool silent = false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_parameters(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_tau( const double new_tau, bool silent = false ) // Truncation parameter
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_tau(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_c( const double new_c, bool silent = false ) // Concentration
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_c(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	const int override_rhmvir( const double &new_rhmvir, bool silent =
			false )
	{
		_rhmvir_cache_ = new_rhmvir;
		hmvir_cached = true;
		return 0;
	}
	const int override_rhmtot( const double &new_rhmtot, bool silent =
			false )
	{
		_rhmtot_cache_ = new_rhmtot;
		hmtot_cached = true;
		return 0;
	}

#endif

#if (1) // Basic get functions
	// All return values in default unit set (m, s, kg, K, rad, C)

	virtual const double mtot() const // Total mass
	{
		return 0;
	}

	virtual const double rt( const bool silent = false ) const // Tidal/truncation radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameters(...) must be overridden to be used.\n";
		return 0;
	}

#endif // end basic get functions

	virtual const unsigned int num_parameters() const
	{
		throw std::runtime_error("ERROR: num_parameters must be overridden if it's used.");
		return 0;
	}

	virtual const int get_parameters( std::vector< double > & parameters,
			const bool silent = false ) const // Returns a set of parameters for this halo (all the variables needed to define it - must be defined for each child)
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameters(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

	virtual const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = false ) const // Returns a set of names of this halo's parameters
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::get_parameter_names(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

#endif // end Virtual functions which must be overwritten if they're going to be used

#if (1) // Virtual functions which should be overwritten if at all possible and if they'll be used
	virtual const double enc_mass( const double &r, const bool silent =
			true ) const; // Mass enclosed with sphere of radius r
#endif

#if (1) // Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

	virtual const double rvir() const // Virial radius (exact definition can be chosen per profile if preferred)
	{
		double virial_density_factor = 200;

		return safe_pow(
				2 * mvir() * Gc / ( pow( H(), 2 ) * virial_density_factor ),
				1. / 3. );
	}
	const double hmvir() const // Half virial mass
	{
		return mvir() / 2;
	}
	const double hmtot() const // Half total mass
	{
		return mtot() / 2;
	}
	virtual const double hm() const // Half mass (which depends on hm_type value)
	{
		return ( _hm_type_ == 0 ? hmvir() : hmtot() );
	}
	virtual const double enc_dens( const double &r,
			const bool silent = false ) const // Mean density enclosed with sphere of radius r
	{
		double r_to_use = max( std::fabs( r ), SMALL_FACTOR );
		return enc_mass( r_to_use, silent )
				/ ( 4. / 3. * pi * std::pow( std::fabs( r_to_use ), 3 ) );
	}
	virtual const double rhmvir( const bool silent = false ) const; // Half-virial-mass radius
	virtual const double rhmtot( const bool silent = false ) const; // Half-total-mass radius
	virtual const double rhm( const bool silent = false ) const // Half-mass radius
	{
		return ( _hm_type_ == 0 ? rhmvir( silent ) : rhmtot( silent ) );
	}
	virtual const double vvir() const // Orbital velocity at rvir
	{
		return _v_from_r( rvir() );
	}
	virtual const double vhmvir( const bool silent = false ) const // Orbital velocity at rhmvir
	{
		return _v_from_r( rhmvir( silent ) );
	}
	virtual const double vhmtot( const bool silent = false ) const // Orbital velocity at rhmtot
	{
		return _v_from_r( rhmvir( silent ) );
	}
	virtual const double vhm( const bool silent = false ) const // Orbital velocity at rhm
	{
		return ( _hm_type_ == 0 ? vhmvir( silent ) : vhmtot( silent ) );
	}
	virtual const double vt( const bool silent = false ) const // Orbital velocity at rt
	{
		return _v_from_r( rt( silent ) );
	}
	const double otvir() const // Orbital period at rvir
	{
		return 2 * pi * rvir() / vvir();
	}
	const double othmvir( const bool silent = false ) const // Orbital period at rhmvir
	{
		return 2 * pi * rhmvir( silent ) / safe_d(vhmvir( silent ));
	}
	const double othmtot( const bool silent = false ) const // Orbital period at rhmtot
	{
		return 2 * pi * rhmtot( silent ) / safe_d(vhmtot( silent ));
	}
	const double othm( const bool silent = false ) const // Orbital period at rhm
	{
		return 2 * pi * rhm( silent ) / safe_d(vhm( silent ));
	}
	const double ott( const bool silent = false ) const // Orbital period at rt
	{
		return 2 * pi * rt( silent ) / safe_d(vt( silent ));
	}

	const int set_hm_type( int new_hm_type ) // Whether default half-mass is half virial, half total, or something else
	{
		_hm_type_ = new_hm_type;
		return 0;
	}

#if (1) // advanced get functions

	const double accel( const double &r,
			const bool silent = false ) const // Gravitational acceleration at radius r
	{
		return r == 0 ? 0 : -Gc * enc_mass( r, silent ) / std::pow( r, 2 );
	}
	virtual const double Daccel( const double &r, const bool silent =
			false ) const; // Derivative of acceleration at radius r

#endif // end advanced get functions

#endif // end Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

#if (1) // Other operations

	virtual const int truncate_to_fraction( const double fraction,
			bool silent = false ) // Adjusts parameters of this class to decrease the mass to fraction of its previous mass. Must be defined for each child
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::truncate_to_fraction(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}

#endif

#endif
};

class tNFW_profile: public density_profile
{
	/**********************************
	 tNFW_profile class
	 -------------

	 A virtual class for an object whose
	 mass fits a truncated NFW profile.
	 See ... .

	 Defined by four parameters:
	 mvir0, z, c, and tau

	 Parent class:
	 density profile

	 Derived classes:
	 tNFW_galaxy

	 **********************************/
private:
#if (1) // Private member variables specific to this class
	double _mvir0_;
	double _c_, _tau_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	tNFW_profile();

	tNFW_profile( const double &init_mvir0, const double init_z,
			const double init_c = -1, const double init_tau = -1 );

#endif // End constructors

	// Destructor
	~tNFW_profile();

#if (1) // Set functions

	const int set_mvir( const double &new_halo_mass, const bool silent =
			false );
	const int set_parameters( const unsigned int num_parameters,
			const std::vector< double > & new_parameters,
			const bool silent = false );

	const int set_z( const double new_z );
	const int set_tau( const double new_halo_tau, const bool silent = false );
	const int set_c( const double new_halo_c, const bool silent = false );

#endif // End set functions

#if (1) //Basic get functions

	const double mvir0() const;

	const double hmvir() const;
	const double mvir() const;
	const double mtot() const;

	const double rvir() const;
	const double rt( const bool silent = false ) const;
	const double rs() const;

	const double vvir() const;

	const double c() const;
	const double tau() const;
#endif // end basic get functions

#if (1) // advanced get functions

	const double dens( const double &r ) const;
	const double enc_mass( const double &r,
			const bool silent = false ) const;
	const unsigned int num_parameters() const
	{
		return 4; // Mass, redshift, c, and tau
	}
	const int get_parameters( std::vector< double > & parameters,
			const bool silent = true ) const;

	const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = true ) const;
#endif // end advanced get functions

#if (1) // Other operations

	virtual const int truncate_to_fraction( const double fraction,
			const bool silent = false );
	virtual density_profile *density_profile_clone() const
	{
		return new tNFW_profile( *this );
	}
	virtual tNFW_profile *tNFW_profile_clone() const
	{
		return new tNFW_profile( *this );
	}

#endif

#endif // end public variables and functions
};

class point_mass_profile: public density_profile
{
	/**********************************
	 point_mass_profile class
	 ------------------------

	 A virtual class for a point mass

	 Defined by two parameters:
	 mass and z

	 Parent class:
	 density profile

	 Derived classes:
	 (none)

	 **********************************/
private:
#if (1) // private member variables specific to this class

	double _mass_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	point_mass_profile();

	point_mass_profile( const double init_mass, const double init_z );

#endif // End constructors

	// Destructor
	~point_mass_profile();

#if (1) // Set functions
	virtual const int set_mvir( const double &new_halo_mass, bool silent =
			false );
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< double > &new_parameters,
			bool silent = false );
#endif // End set functions

#if (1) // Basic get functions
	const double mass() const;

	const double hmvir() const;
	const double mvir() const;
	const double mtot() const;

	const double rvir() const;
	const double rt(const bool silent = false) const;
	const double rs() const;

	const double vvir() const;

#endif // end basic get functions

#if (1) // advanced get functions
	const double dens( const double &r ) const;
	const double enc_dens( const double &r,
			const bool silent = false ) const;
	const double enc_mass( const double &r, const bool silent =
				true ) const; // Mass enclosed with sphere of radius r
	const unsigned int num_parameters() const
	{
		return 2; // Mass and redshift
	}
	const int get_parameters( std::vector< double > & parameters,
			const bool silent = false ) const;

	const int get_parameter_names( std::vector< std::string > & parameter_names,
			const bool silent = false ) const;
#endif // end advanced get functions

#if (1) // Other operations

	const int truncate_to_fraction( const double fraction,
			bool silent = false );
	virtual density_profile *density_profile_clone() const
	{
		return new point_mass_profile( *this );
	}
	virtual point_mass_profile *point_mass_profile_clone() const
	{
		return new point_mass_profile( *this );
	}

#endif

#endif // end public variables and functions
};
// class point_mass_profile

class accel_function: public functor< double >
{
	/**********************************
	 accel_function class
	 -----------------------------

	 Function class for acceleration within a density profile

	 Parent class: function_class (from brg_functions)

	 **********************************/
private:
	const density_profile *_host_ptr_;
public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;

	accel_function();
	accel_function( const density_profile *init_host_ptr );
	virtual ~accel_function()
	{
	}
};
// class accel_function

class solve_rhm_function: public functor< double >
{
	/**********************************
	 solve_rhm_function class
	 -----------------------------

	 Function class for solving the half-mass
	 radius of a halo.

	 Parent class: function_class (from brg_functions)

	 **********************************/
private:
	const density_profile *_host_ptr_;double _target_mass_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}
	const int set_target_mass( const double &new_target_mass );
	const double & target_mass()
	{
		return _target_mass_;
	}

	const int operator ()( const double & in_param,
	double & out_param, const bool silent = false ) const;

	solve_rhm_function();
	solve_rhm_function( const density_profile *init_host,
			const double &init_target_mass );

};
// end class unitless_solve_rhm_function

class spherical_density_function: public functor< double >
{
	/**********************************
	 spherical_density_function class
	 -----------------------------

	 Function class integrating density in a sphere

	 Parent class: function_class (from brg_functions)

	 **********************************/
private:
	const density_profile *_host_ptr_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int operator()( const double & in_param,
	double & out_param, const bool silent = false ) const;

	spherical_density_function();
	spherical_density_function( const density_profile *init_host );
	virtual ~spherical_density_function()
	{
	}
};

#endif // end Class Definitions

} // end namespace SALTSA

#endif

