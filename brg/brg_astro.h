/**********************************************************************\
brg_astro.h
 -----------

 If this header is used, the source file brg_astro.cpp must be included
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
 in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_ASTRO_H_INCLUDED__
#define __BRG_ASTRO_H_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include "brg_units.h"
#include "brg_functions.h"
#include "brg_cache.hpp"

/** Constant Definitions **/
#if (1)
/**********************************************************************
 This section stores the definitions of various physical and
 astrophysical constants, stored as variables. Key points:

 -H_0 is set to 70 km/s/Mpc. If this is not changed, all results should
 be assumed to be in h_70 units.
 -Cosmological values are based on best-fit WMAP + priors values
 from Hinshaw et al. 2012.
 -If the _BRG_USE_UNITS_ parameter is set (see the brg_global.h file), values
 that have units will be declared as unit_objs.

 \**********************************************************************/

// Physical Constants
#if (1)
#ifndef __PI_DEFINED__
#define __PI_DEFINED__
const double pi = 3.14159265358979323846;
#endif

#ifndef __PHYS_VARS_DEFINED__
#define __PHYS_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj Gc(6.67384e-11,3,-2,-1,0,0); // In m^3 s^-2 kg^-1
#else
const float Gc = 6.67384e-11; // In m^3 s^-2 kg^-1
#endif
const BRG_VELOCITY c = unitconv::ctomps;

#endif

#endif // end Physical Constants

// Cosmological Parameters
#if (1)
#ifndef __COSMO_VARS_DEFINED__
#define __COSMO_VARS_DEFINED__

#ifdef _BRG_USE_UNITS_
const brgastro::unit_obj H_0(70*unitconv::kmtom/unitconv::stos/unitconv::Mpctom,0,-1,0,0,0,0); // So all results will implicitly be in h_70 units
#else
const double H_0 = 70 * unitconv::kmtom / unitconv::stos / unitconv::Mpctom; // So all results will implicitly be in h_70 units
#endif

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

namespace brgastro
{

/** Class Forward Declarations **/
#if (1)
/**********************************************************************
 This section forward declares all class defined in this header, allowing
 them to point to each other if necessary.

 \**********************************************************************/

class redshift_obj;
// Concrete base

class sky_obj;
// Abstract base
class galaxy;
// Concrete
class group;
// Concrete

class density_profile;
// Abstract base
class tNFW_profile;
// Concrete
class point_mass_profile;
// Concrete

class tNFW_galaxy;
// Concrete
class tNFW_group;
// Concrete

class spherical_density_function;
// Concrete
class projected_density_function;
// Concrete
class cylindrical_density_function;
// Concrete
class accel_function;
// Concrete
class solve_rhm_function;
// Concrete

class offset_ring_dens_function;
// Concrete
class offset_circ_dens_function;
// Concrete
class offset_WLsig_function;
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

class grid_cache
{
private:
	static int _ra_grid_change_num_, _dec_grid_change_num_,
			_z_grid_change_num_;
	static BRG_ANGLE _ra_grid_min_, _ra_grid_max_, _ra_grid_step_;
	static BRG_ANGLE _dec_grid_min_, dec_grid_max_val, _dec_grid_step_;
	static double _z_grid_min_, _z_grid_max_, _z_grid_step_;
public:
	// Set functions
#if (1)
	const int set_ra_grid( const BRG_ANGLE new_ra_grid_min,
			const BRG_ANGLE new_ra_grid_max, const BRG_ANGLE new_ra_grid_step )
	{
		_ra_grid_min_ = new_ra_grid_min;
		_ra_grid_max_ = new_ra_grid_max;
		_ra_grid_step_ = new_ra_grid_step;
		_ra_grid_change_num_++;
		return 0;
	}

	const int set_dec_grid( const BRG_ANGLE new_dec_grid_min,
			const BRG_ANGLE new_dec_grid_max,
			const BRG_ANGLE new_dec_grid_step )
	{
		_dec_grid_min_ = new_dec_grid_min;
		dec_grid_max_val = new_dec_grid_max;
		_dec_grid_step_ = new_dec_grid_step;
		_dec_grid_change_num_++;
		return 0;
	}

	const int set_z_grid( const double new_z_grid_min,
			const double new_z_grid_max, const double new_z_grid_step )
	{
		_z_grid_min_ = new_z_grid_min;
		_z_grid_max_ = new_z_grid_max;
		_z_grid_step_ = new_z_grid_step;
		_z_grid_change_num_++;
		return 0;
	}
#endif

	// Get functions
#if (1)

	const int ra_grid_change_num()
	{
		return _ra_grid_change_num_;
	}
	const int dec_grid_change_num()
	{
		return _dec_grid_change_num_;
	}
	const int z_grid_change_num()
	{
		return _z_grid_change_num_;
	}
	const BRG_ANGLE ra_grid_min()
	{
		return _ra_grid_min_;
	}
	const BRG_ANGLE dec_grid_min()
	{
		return _dec_grid_min_;
	}
	const double z_grid_min()
	{
		return _z_grid_min_;
	}
	const BRG_ANGLE ra_grid_max()
	{
		return _ra_grid_max_;
	}
	const BRG_ANGLE dec_grid_max()
	{
		return dec_grid_max_val;
	}
	const double z_grid_max()
	{
		return _z_grid_max_;
	}
	const BRG_ANGLE ra_grid_step()
	{
		return _ra_grid_step_;
	}
	const BRG_ANGLE dec_grid_step()
	{
		return _dec_grid_step_;
	}
	const double z_grid_step()
	{
		return _z_grid_step_;
	}
#endif

};
// class grid_cache

class dfa_cache : public brg_cache<dfa_cache>
{
	// "Distance from angle" cache
private:

	DECLARE_BRG_CACHE_STATIC_VARS();

	friend class brg_cache;

protected:

	const std::string _name_base() const throw()
	{
		return "dfa";
	}

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class tfa_cache

class add_cache
{
	// "Angular Diameter Distance" cache
	static bool loaded;
	static std::string file_name;
	static double z_min, z_max, z_step;
	static int z_res;
	static std::vector< std::vector< double > > dmod;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_min, const double new_z_max,
			const double new_z_step, const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );
	const BRG_UNITS get( const double z1, const double z2, const bool silent =
			false );

};
// class add_cache

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

#ifdef _BRG_USE_UNITS_

	// Tells what units the result should have. Only the units matter in the return, not the value
	const brgastro::unit_obj _units() const throw()
	{
		return brgastro::unit_obj(0,0,1,0,0,0,0);
	}
	const brgastro::unit_obj _inverse_units() const throw()
	{
		return brgastro::unit_obj(0);
	}

#endif // _BRG_USE_UNITS_

	// Long-form calculation function.
	const int _calculate( const double in_params, double & out_params ) const;

public:

};
// class tfa_cache

class tNFW_sig_cache
{
	// Weak lensing from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static std::vector< std::vector< std::vector< double > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r_halo, const bool silent = false ); // Takes in standard units

};
// class tNFW_sig_cache

class tNFW_offset_sig_cache
{
	// Offset weak lensing signal from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static double offset_r_min, offset_r_max, offset_r_step;
	static int offset_r_res;
	static std::vector< std::vector< std::vector< std::vector< double > > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const double new_offset_r_min, const double new_offset_r_max,
			const double new_offset_r_step, const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r, const BRG_DISTANCE offset_r,
			const bool silent = false ); // Takes in standard units

};
// class tNFW_offset_sig_cache

class tNFW_group_sig_cache
{
	// Group weak lensing signal from tNFW profile cache
	static bool loaded;
	static std::string file_name;
	static double halo_z_min, halo_z_max, halo_z_step;
	static int halo_z_res;
	static double r_min, r_max, r_step;
	static int r_res;
	static double halo_m_min, halo_m_max, halo_m_step;
	static int halo_m_res;
	static double group_c_min, group_c_max, group_c_step;
	static int group_c_res;
	static std::vector< std::vector< std::vector< std::vector< double > > > > signal;
	static std::string header_string;
	static int sig_digits;
	const int load( const bool silent = false );
	const int unload();
	const int calc( const bool silent = false );
	const int output( const bool silent = false );
public:
	const int set_file_name( const std::string new_name );
	const int set_range( const double new_z_halo_min,
			const double new_z_halo_max, const double new_z_halo_step,
			const double new_m_halo_min, const double new_m_halo_max,
			const double new_m_halo_step, const double new_r_min,
			const double new_r_max, const double new_r_step,
			const double new_group_c_min, const double new_group_c_max,
			const double new_group_c_step, const bool silent = false );
	const int set_precision( const int new_precision,
			const bool silent = false );

	const BRG_UNITS get( const double z_halo, const BRG_MASS m_halo,
			const BRG_DISTANCE r, const double group_c, const bool silent =
					false ); // Takes in standard units

};
// class tNFW_group_sig_cache

#endif // End static class definitions

/** Function Declarations **/
#if (1)
/**********************************************************************
 This section defines various functions which I find useful for
 astrophysical calculations. All are declared in the namespace brgastro.

 \**********************************************************************/

// Functions to get grid integers or grid boundaries from integers
const int get_ra_grid( const BRG_ANGLE &ra );
const int get_dec_grid( const BRG_ANGLE &dec );
const int get_z_grid( const double z );

const BRG_ANGLE get_ra_grid_lower( const int ra_grid );
const BRG_ANGLE get_dec_grid_lower( const int dec_grid );
const double get_z_grid_lower( const int z_grid );

const BRG_ANGLE get_ra_grid_upper( const int ra_grid );
const BRG_ANGLE get_dec_grid_upper( const int dec_grid );
const double get_z_grid_upper( const int z_grid );

const BRG_ANGLE get_ra_grid_mid( const int ra_grid );
const BRG_ANGLE get_dec_grid_mid( const int dec_grid );
const double get_z_grid_mid( const int z_grid );

// Functions to get transverse distance (in m) from angle (in rad) or vice-versa
const BRG_DISTANCE dfa( const BRG_ANGLE &da, const double z );
const BRG_DISTANCE dfa( const BRG_ANGLE &a1, const BRG_ANGLE &a2,
		const double z );
const BRG_DISTANCE dfa( const BRG_ANGLE &a1x, const BRG_ANGLE &a1y,
		const BRG_ANGLE &a2x, const BRG_ANGLE &a2y, const double z );
const BRG_DISTANCE dfa( const sky_obj *obj1, const sky_obj *obj2,
		const double z = -1 );

const BRG_ANGLE afd( const BRG_DISTANCE &dd, const double z );
const BRG_ANGLE afd( const BRG_DISTANCE &d1, const BRG_DISTANCE &d2,
		const double z );
const BRG_ANGLE afd( const BRG_DISTANCE &d1x, const BRG_DISTANCE &d1y,
		const BRG_DISTANCE &d2x, const BRG_DISTANCE &d2y, const double z );

// Functions to work between redshift, scale factor, and time (in s, with zero = present day)
const double zfa( const double a );
const double afz( const double z );

const BRG_TIME tfz( const double z );
const BRG_TIME tfa( const double z );
const double zft( const BRG_TIME &t );
const double aft( const BRG_TIME &t );

// Functions to integrate out distances
const double integrate_add( const double z1, const double z2 );
const double integrate_cmd( const double z1, const double z2 );
const double integrate_Ld( const double z1, const double z2 );
const double integrate_ltd( const double z1, const double z2 );
const double integrate_add( const double z );
const double integrate_cmd( const double z );
const double integrate_Ld( const double z );
const double integrate_ltd( const double z );
const double integrate_distance( const double z1, const double z2,
		const int mode, const int resolution = 10000 );

// Functions relating to tNFW profiles
inline const double cfm( const BRG_MASS mass, const double z = 0 ) // Concentration from mass relationship, from Neto
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
const double taufm( const double mratio, double c, double tau_init = 0, double precision=0.00001,
		const bool silent = false ); // tau from Mtot/Mvir
inline const double delta_c( const double conc ) // Simple function of concentration used as a step in calculating NFW densities
{
	return ( 200. / 3. ) * std::pow( conc, 3 )
			/ ( log( 1 + conc ) - conc / ( 1 + conc ) );
}

// Function to estimate orbital period from current position and velocity in a density profile
// Note that this is merely an estimate from analogy to calculations in a Keplerian potential
const BRG_TIME period( const density_profile *host, const BRG_DISTANCE &r,
		const BRG_VELOCITY &vr, const BRG_VELOCITY &vt = 0 );

// Lensing functions
const BRG_DISTANCE ad_distance( double z1, double z2 = 0 );
const BRG_UNITS sigma_crit( const double z_lens, const double z_source );

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
	mutable int _z_grid_;
	mutable bool _z_grid_cached_;
	mutable int _local_z_grid_change_num_;
public:
	// Constructor
	redshift_obj( const double init_z = 0, const double init_z_err = 0 )
	{
		_z_ = init_z;
		_z_err_ = init_z_err;
		_z_grid_ = 0;
		_z_grid_cached_ = false;
		_local_z_grid_change_num_ = -1;
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
		_z_grid_cached_ = false;
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
	virtual const int z_grid() const;
#endif

	// H(z) at either the redshift of the object or a specified redshift, given in units of m/s^2
	virtual const BRG_UNITS H( double ztest = -1 ) const;

	// Clone function - Not needed in current implementation
	// virtual redshift_obj *redshift_obj_clone()=0;
};

class sky_obj: public virtual redshift_obj
{
	/**********************************
	 sky_obj class
	 -------------

	 An abstract class for anything in the sky.

	 Derived classes:
	 galaxy
	 group

	 **********************************/
private:
#if (1)
	std::string _ID_; // Name for it or ID number

	int _index_; // Position in an array

	double _weight_;

	mutable int _ra_grid_, _dec_grid_;

	mutable bool _ra_grid_cached_, _dec_grid_cached_;
	mutable int _local_ra_grid_change_num_, _local_dec_grid_change_num_;

	BRG_ANGLE _ra_, _ra_err_;BRG_ANGLE _dec_, _dec_err_;
#endif
public:
#if (1)
	// Public member functions

	// Constructor
	sky_obj( BRG_ANGLE init_ra = 0, BRG_ANGLE init_dec = 0, double init_z = 0,
	BRG_ANGLE init_ra_err = 0, BRG_ANGLE init_dec_err = 0, double init_z_err =
			0 ); // Normal constructor

	//sky_obj(sky_obj other_sky_obj); // Copy constructor (Default is fine for us)

	virtual ~sky_obj()
	{
	} // Virtual destructor, in case it needs to be overridden

	virtual const int clear(); // Resets all variables to zero
	virtual const int partial_clear(); // Resets all variables which can't be initialized

#if (1) // Set functions
	virtual const int set_ra( const BRG_ANGLE &new_ra );
	virtual const int set_dec( const BRG_ANGLE &new_dec );
	virtual const int set_ra_err( const BRG_ANGLE &new_ra_err );
	virtual const int set_dec_err( const BRG_ANGLE &new_dec_err );
	virtual const int set_ra_dec( const BRG_ANGLE &new_ra,
			const BRG_ANGLE &new_dec ); // Sets ra and dec
	virtual const int set_ra_dec_z( const BRG_ANGLE &new_ra,
			const BRG_ANGLE &new_dec, const double new_z ); // Sets all values
	virtual const int set_ra_dec_z_err( const BRG_ANGLE &new_ra,
			const BRG_ANGLE &new_dec, const double new_z,
			const BRG_ANGLE &new_ra_err, const BRG_ANGLE &new_dec_err,
			const double new_z_err ); // Sets all values and error
	virtual const int set_ra_dec_err( const BRG_ANGLE &new_ra,
			const BRG_ANGLE &new_dec, const BRG_ANGLE &new_ra_err,
			const BRG_ANGLE &new_dec_err ); // Sets ra and dec and error

	virtual const int set_weight( const double new_weight );
	virtual const int set_index( const int new_index );
	virtual const int set_ID( const std::string &new_ID );
#endif // end set functions

#if (1) //Get functions

	virtual const BRG_ANGLE ra() const
	{
		return _ra_;
	}
	virtual const BRG_ANGLE dec() const
	{
		return _dec_;
	}
	virtual const BRG_ANGLE ra_err() const
	{
		return _ra_err_;
	}
	virtual const BRG_ANGLE dec_err() const
	{
		return _dec_err_;
	}

	virtual const double weight() const
	{
		return _weight_;
	}
	virtual const int index() const
	{
		return _index_;
	}
	virtual const std::string ID() const
	{
		return _ID_;
	}

	virtual const int ra_grid() const;
	virtual const int dec_grid() const;

#endif // end get functions

	// Clone function (to enable copies to be made of pointed-to objects)
	virtual redshift_obj *redshift_obj_clone()=0;
	virtual sky_obj *sky_obj_clone()=0;

#endif

};
// class sky_obj

class galaxy: public sky_obj
{
	/**********************************
	 galaxy class
	 -------------

	 A class for galaxies.

	 Parent class: sky_obj

	 Derived classes: tNFW_galaxy

	 **********************************/

private:
	// private member variables (none needed in current implementation. Should make all private for consistency in future though)

public:

	// Public member variables

	BRG_MASS stellar_mass;
	double umag, gmag, rmag, imag, zmag;
	double umag_err, gmag_err, rmag_err, imag_err, zmag_err;

	double z_phot, z_phot_err;
	double odds;
	double phot_template;

	group *host_group; // Pointer to group this galaxy resides within
	int host_group_index;

	// Public member functions

	// Constructor
	galaxy();

	// Copy constructor
	//galaxy(const galaxy other_galaxy); // Implicit is fine for us, though watch in case we need to make it virtual

	// Virtual destructor
	virtual ~galaxy()
	{
	}
	;

	// Clear function
	virtual const int clear();

	// Clone functions
	virtual redshift_obj *redshift_obj_clone()
	{
		return new galaxy( *this );
	}
	virtual sky_obj *sky_obj_clone()
	{
		return new galaxy( *this );
	}
	virtual galaxy *galaxy_clone()
	{
		return new galaxy( *this );
	}

};
// class unitless_galaxy

class group: public sky_obj
{
	/**********************************
	 group class
	 -------------

	 A class for groups and clusters of galaxies.

	 Parent class: sky_obj

	 Derived classes:
	 (none)

	 **********************************/

private:

public:
	// Public member variables

	double z_phot, z_phot_err;
	double odds;

	int BCG_index;

	int num_members;
	std::vector< int > member_indices;
	std::vector< galaxy * > members;

	// Public member functions

	// Constructor
	group( int init_num_members = 0 );
	group( double init_mass, double init_z, double init_c = -1,
			double init_tau = -1 );

	// Copy constructor
	group( const group &unitless_group );

	// Destructor
	virtual ~group();

	// Functions to set parameters
	virtual const int clear();
	virtual const int resize( int new_num_members );
	virtual const int set_member( int index, galaxy * new_member,
			const bool silent = false );
	virtual const int set_member_index( int index, int new_member_index,
			const bool silent = false );
	virtual const int add_member( galaxy * new_member, const bool silent =
			false );
	virtual const int remove_member( galaxy * rem_member, const bool silent =
			false );

	// Clone functions
	virtual redshift_obj *redshift_obj_clone()
	{
		return new group( *this );
	}
	virtual sky_obj *sky_obj_clone()
	{
		return new group( *this );
	}
	virtual group *group_clone()
	{
		return new group( *this );
	}

};
// class group

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

	mutable BRG_DISTANCE _rhmvir_cache_, _rhmtot_cache_;

	BRG_VELOCITY _v_from_r( BRG_DISTANCE r ) const
	{
		BRG_UNITS a;

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

	virtual const BRG_MASS mvir() const =0; // Virial mass (exact definition can be chosen per profile)
	virtual const BRG_UNITS dens( const BRG_DISTANCE &r ) const =0; // Local density at radius r

	virtual density_profile *density_profile_clone() const =0; // Creates a clone of this
#endif

#if (1) // Virtual functions which must be overwritten if they're going to be used

#if (1) // Set functions - will return 1 if profile doesn't support this method of setting
	// All take values in default unit set (m, s, kg, K, rad, C)
	virtual const int set_mvir( const BRG_MASS &new_mvir, bool silent = false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_mvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_vvir( const BRG_VELOCITY &new_vvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_vvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rvir( const BRG_DISTANCE &new_rvir, bool silent =
			false )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rvir(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rs( const BRG_DISTANCE &new_rs, bool silent = false ) // Scale radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rs(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_rt( const BRG_DISTANCE &new_rt, bool silent = false ) // Tidal/truncation radius
	{
		if ( !silent )
			std::cerr
					<< "ERROR: density_profile::set_rt(...) must be overridden to be used.\n";
		return MUST_OVERRIDE_ERROR;
	}
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > &parameters, bool silent = false )
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
	const int override_rhmvir( const BRG_DISTANCE &new_rhmvir, bool silent =
			false )
	{
		_rhmvir_cache_ = new_rhmvir;
		hmvir_cached = true;
		return 0;
	}
	const int override_rhmtot( const BRG_DISTANCE &new_rhmtot, bool silent =
			false )
	{
		_rhmtot_cache_ = new_rhmtot;
		hmtot_cached = true;
		return 0;
	}

#endif

#if (1) // Basic get functions
	// All return values in default unit set (m, s, kg, K, rad, C)

	virtual const BRG_MASS mtot() const // Total mass
	{
		return 0;
	}

	virtual const BRG_DISTANCE rt( const bool silent = false ) const // Tidal/truncation radius
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

	virtual const int get_parameters( std::vector< BRG_UNITS > & parameters,
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
	virtual const BRG_MASS enc_mass( const BRG_DISTANCE &r, const bool silent =
			true ) const; // Mass enclosed with sphere of radius r
	virtual const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Projected surface density at radius R
	virtual const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R,
			const bool silent = true ) const; // Mass enclosed within a cylinder of radius R
	virtual const BRG_UNITS offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from position offset by offset_R
	virtual const BRG_UNITS group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const; // Expected weak lensing signal in tangential shear Delta-Sigma at radius R from ensemble of satellites in group with satellite concentration group_c
	virtual const BRG_UNITS quick_WLsig( const BRG_DISTANCE &R,
			const bool silent = true ) const // As deltasigma, but uses cache to speed it up if overwritten
	{
		return deltasigma( R, silent );
	}
	virtual const BRG_UNITS quick_offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = true ) const // As offset_WLsig, but uses cache to speed it up if overwritten
	{
		return offset_WLsig( R, offset_R, silent );
	}
	virtual const BRG_UNITS semiquick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const // As group_WLsig, but uses offset_WLsig cache to speed it up if overwritten
	{
		return group_WLsig( R, group_c, silent );
	}
	virtual const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = true ) const // As deltasigma, but uses group_WLsig cache to speed it up if overwritten
	{
		return group_WLsig( R, group_c, silent );
	}
#endif

#if (1) // Virtual functions which shouldn't be overwritten in most cases and non-virtual functions

	virtual const BRG_DISTANCE rvir() const // Virial radius (exact definition can be chosen per profile if preferred)
	{
		double virial_density_factor = 200;

		return safe_pow(
				2 * mvir() * Gc / ( pow( H(), 2 ) * virial_density_factor ),
				1. / 3. );
	}
	const BRG_MASS hmvir() const // Half virial mass
	{
		return mvir() / 2;
	}
	const BRG_MASS hmtot() const // Half total mass
	{
		return mtot() / 2;
	}
	virtual const BRG_MASS hm() const // Half mass (which depends on hm_type value)
	{
		return ( _hm_type_ == 0 ? hmvir() : hmtot() );
	}
	virtual const BRG_UNITS enc_dens( const BRG_DISTANCE &r,
			const bool silent = false ) const // Mean density enclosed with sphere of radius r
	{
		BRG_DISTANCE r_to_use = max( std::fabs( r ), SMALL_FACTOR );
		return enc_mass( r_to_use, silent )
				/ ( 4. / 3. * pi * std::pow( std::fabs( r_to_use ), 3 ) );
	}
	virtual const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const // Mean surface density enclosed within a cylinder of radius R
	{
		BRG_DISTANCE R_to_use = max( std::fabs( R ), SMALL_FACTOR );
		return proj_enc_mass( R_to_use, silent )
				/ ( pi * std::pow( std::fabs( R_to_use ), 2 ) );
	}
	virtual const BRG_DISTANCE rhmvir( const bool silent = false ) const; // Half-virial-mass radius
	virtual const BRG_DISTANCE rhmtot( const bool silent = false ) const; // Half-total-mass radius
	virtual const BRG_DISTANCE rhm( const bool silent = false ) const // Half-mass radius
	{
		return ( _hm_type_ == 0 ? rhmvir( silent ) : rhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vvir() const // Orbital velocity at rvir
	{
		return _v_from_r( rvir() );
	}
	virtual const BRG_VELOCITY vhmvir( const bool silent = false ) const // Orbital velocity at rhmvir
	{
		return _v_from_r( rhmvir( silent ) );
	}
	virtual const BRG_VELOCITY vhmtot( const bool silent = false ) const // Orbital velocity at rhmtot
	{
		return _v_from_r( rhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vhm( const bool silent = false ) const // Orbital velocity at rhm
	{
		return ( _hm_type_ == 0 ? vhmvir( silent ) : vhmtot( silent ) );
	}
	virtual const BRG_VELOCITY vt( const bool silent = false ) const // Orbital velocity at rt
	{
		return _v_from_r( rt( silent ) );
	}
	const BRG_TIME otvir() const // Orbital period at rvir
	{
		return 2 * pi * rvir() / vvir();
	}
	const BRG_TIME othmvir( const bool silent = false ) const // Orbital period at rhmvir
	{
		return 2 * pi * rhmvir( silent ) / safe_d(vhmvir( silent ));
	}
	const BRG_TIME othmtot( const bool silent = false ) const // Orbital period at rhmtot
	{
		return 2 * pi * rhmtot( silent ) / safe_d(vhmtot( silent ));
	}
	const BRG_TIME othm( const bool silent = false ) const // Orbital period at rhm
	{
		return 2 * pi * rhm( silent ) / safe_d(vhm( silent ));
	}
	const BRG_TIME ott( const bool silent = false ) const // Orbital period at rt
	{
		return 2 * pi * rt( silent ) / safe_d(vt( silent ));
	}

	const int set_hm_type( int new_hm_type ) // Whether default half-mass is half virial, half total, or something else
	{
		_hm_type_ = new_hm_type;
		return 0;
	}

#if (1) // advanced get functions

	const BRG_UNITS accel( const BRG_DISTANCE &r,
			const bool silent = false ) const // Gravitational acceleration at radius r
	{
		return r == 0 ? 0 : -Gc * enc_mass( r, silent ) / std::pow( r, 2 );
	}
	virtual const BRG_UNITS Daccel( const BRG_DISTANCE &r, const bool silent =
			false ) const; // Derivative of acceleration at radius r
	const BRG_UNITS deltasigma( const BRG_DISTANCE &R, const bool silent =
			false ) const // Weak lensing signal in tangential shear Delta-Sigma at radius R
	{
		return proj_enc_dens( R, silent ) - proj_dens( R, silent );
	}

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
	BRG_MASS _mvir0_;
	double _c_, _tau_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	tNFW_profile();

	tNFW_profile( const BRG_MASS &init_mvir0, const double init_z,
			const double init_c = -1, const double init_tau = -1 );

#endif // End constructors

	// Destructor
	~tNFW_profile();

#if (1) // Set functions

	const int set_mvir( const BRG_MASS &new_halo_mass, const bool silent =
			false );
	const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > & new_parameters,
			const bool silent = false );

	const int set_z( const double new_z );
	const int set_tau( const double new_halo_tau, const bool silent = false );
	const int set_c( const double new_halo_c, const bool silent = false );

#endif // End set functions

#if (1) //Basic get functions

	const BRG_MASS mvir0() const;

	const BRG_MASS hmvir() const;
	const BRG_MASS mvir() const;
	const BRG_MASS mtot() const;

	const BRG_DISTANCE rvir() const;
	const BRG_DISTANCE rt( const bool silent = false ) const;
	const BRG_DISTANCE rs() const;

	const BRG_VELOCITY vvir() const;

	const double c() const;
	const double tau() const;
#endif // end basic get functions

#if (1) // advanced get functions

	const BRG_UNITS dens( const BRG_DISTANCE &r ) const;
	const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const BRG_MASS enc_mass( const BRG_DISTANCE &r,
			const bool silent = false ) const;
	const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_WLsig( const BRG_DISTANCE &R, const bool silent =
			false ) const;
	const BRG_UNITS quick_offset_WLsig( const BRG_DISTANCE &R,
			const BRG_DISTANCE &offset_R, const bool silent = false ) const;
	const BRG_UNITS semiquick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = false ) const;
	const BRG_UNITS quick_group_WLsig( const BRG_DISTANCE &R,
			const double group_c, const bool silent = false ) const;
	const unsigned int num_parameters() const
	{
		return 4; // Mass, redshift, c, and tau
	}
	const int get_parameters( std::vector< BRG_UNITS > & parameters,
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

	BRG_MASS _mass_;

#endif // end private member variables

public:
#if (1) // public variables and functions

#if (1) // Constructors
	point_mass_profile();

	point_mass_profile( const BRG_MASS init_mass, const double init_z );

#endif // End constructors

	// Destructor
	~point_mass_profile();

#if (1) // Set functions
	virtual const int set_mvir( const BRG_MASS &new_halo_mass, bool silent =
			false );
	virtual const int set_parameters( const unsigned int num_parameters,
			const std::vector< BRG_UNITS > &new_parameters,
			bool silent = false );
#endif // End set functions

#if (1) // Basic get functions
	const BRG_MASS mass() const;

	const BRG_MASS hmvir() const;
	const BRG_MASS mvir() const;
	const BRG_MASS mtot() const;

	const BRG_DISTANCE rvir() const;
	const BRG_DISTANCE rt(const bool silent = false) const;
	const BRG_DISTANCE rs() const;

	const BRG_VELOCITY vvir() const;

#endif // end basic get functions

#if (1) // advanced get functions
	const BRG_UNITS dens( const BRG_DISTANCE &r ) const;
	const BRG_UNITS proj_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const BRG_UNITS enc_dens( const BRG_DISTANCE &r,
			const bool silent = false ) const;
	const BRG_MASS enc_mass( const BRG_DISTANCE &r, const bool silent =
				true ) const; // Mass enclosed with sphere of radius r
	const BRG_UNITS proj_enc_dens( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const BRG_MASS proj_enc_mass( const BRG_DISTANCE &R,
			const bool silent = false ) const;
	const unsigned int num_parameters() const
	{
		return 2; // Mass and redshift
	}
	const int get_parameters( std::vector< BRG_UNITS > & parameters,
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

class tNFW_galaxy: public tNFW_profile, public galaxy
{
	// Simple combination of the two classes. Only complicated part is the z value
public:
	tNFW_galaxy()
	{
		galaxy();
		tNFW_profile();
	}
	virtual density_profile *density_profile_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual tNFW_profile *tNFW_profile_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual galaxy *galaxy_clone() const
	{
		return new tNFW_galaxy( *this );
	}
	virtual sky_obj *sky_obj_clone() const
	{
		return new tNFW_galaxy( *this );
	}
};

class accel_function: public functor< BRG_UNITS >
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

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	accel_function();
	accel_function( const density_profile *init_host_ptr );
	virtual ~accel_function()
	{
	}
};
// class accel_function

class solve_rhm_function: public functor< BRG_UNITS >
{
	/**********************************
	 solve_rhm_function class
	 -----------------------------

	 Function class for solving the half-mass
	 radius of a halo.

	 Parent class: function_class (from brg_functions)

	 **********************************/
private:
	const density_profile *_host_ptr_;BRG_MASS _target_mass_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}
	const int set_target_mass( const BRG_MASS &new_target_mass );
	const BRG_MASS & target_mass()
	{
		return _target_mass_;
	}

	const int operator ()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	solve_rhm_function();
	solve_rhm_function( const density_profile *init_host,
			const BRG_MASS &init_target_mass );

};
// end class unitless_solve_rhm_function

class spherical_density_function: public functor< BRG_UNITS >
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

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	spherical_density_function();
	spherical_density_function( const density_profile *init_host );
	virtual ~spherical_density_function()
	{
	}
};

class projected_density_function: public functor< BRG_UNITS >
{
	/**********************************
	 projected_density_function class
	 -----------------------------

	 Function class integrating density along a projected line

	 Parent class: function_class (from brg_functions)

	 **********************************/
private:
	const density_profile *_host_ptr_;BRG_UNITS _offset_R_;

public:

	const int set_host_ptr( const density_profile *new_host );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_offset_R( const BRG_DISTANCE &new_offset_R );
	const BRG_DISTANCE offset_R()
	{
		return _offset_R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	projected_density_function();
	projected_density_function( const density_profile *init_host,
			const BRG_DISTANCE &init_offset_R );
	virtual ~projected_density_function()
	{
	}
};

class cylindrical_density_function: public functor< BRG_UNITS >
{
	/**********************************
	 cylindrical_density_function class
	 -----------------------------

	 Function class integrating density in a cylinder

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

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	cylindrical_density_function();
	cylindrical_density_function( const density_profile *init_host );
	virtual ~cylindrical_density_function()
	{
	}
};

class offset_ring_dens_function: public functor< BRG_UNITS >
{

	const density_profile *_host_ptr_;BRG_DISTANCE _R0_, _R_;
public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R0( const BRG_DISTANCE &new_R_0 );
	const BRG_DISTANCE & R0()
	{
		return _R0_;
	}

	const int set_R( const BRG_DISTANCE &new_R );
	const BRG_DISTANCE & R()
	{
		return _R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_ring_dens_function();
	offset_ring_dens_function( const density_profile *new_host,
			const BRG_DISTANCE &new_R_0 = 0, const BRG_DISTANCE &new_R = 0 );

};

class offset_circ_dens_function: public functor< std::vector< BRG_UNITS > >
{
private:
	const density_profile *_host_ptr_;BRG_DISTANCE _R0_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R0( const BRG_DISTANCE &new_R0 );
	const BRG_DISTANCE & R0()
	{
		return _R0_;
	}

	const int operator()( const std::vector< BRG_UNITS > & in_params,
			std::vector< BRG_UNITS > & out_params,
			const bool silent = false ) const;

	offset_circ_dens_function();
	offset_circ_dens_function( const density_profile *new_host,
			const BRG_DISTANCE &new_R_0 = 0 );
};

class offset_WLsig_function: public functor< BRG_UNITS >
{

private:

	const density_profile *_host_ptr_;BRG_DISTANCE _R_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_R( const BRG_DISTANCE &new_R );
	const BRG_DISTANCE & R()
	{
		return _R_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_WLsig_function();
	offset_WLsig_function( const density_profile *init_host,
			const BRG_DISTANCE &new_R = 0 );

};

class offset_WLsig_weight_function: public functor< BRG_UNITS >
{

private:

	const density_profile *_host_ptr_;
	double _c_;

public:

	const int set_host_ptr( const density_profile *new_host_ptr );
	const density_profile * host_ptr()
	{
		return _host_ptr_;
	}

	const int set_c( const double new_c );
	const double c()
	{
		return _c_;
	}

	const int operator()( const BRG_UNITS & in_param,
	BRG_UNITS & out_param, const bool silent = false ) const;

	offset_WLsig_weight_function();
	offset_WLsig_weight_function( const density_profile *new_host,
			const double init_c = -1 );

};

#endif // end Class Definitions

} // end namespace brgastro

#endif

