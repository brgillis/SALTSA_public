#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <new>
#include <fstream>
#include "brg_units.h"
#include "brg_functions.h"
#include "brg_astro.h"
#include "brg_calculus.hpp"
#include "brg_solvers.hpp"

using namespace std;

/** Static Class Initialisation **/
#if (1)
// Initialisation for brgastro::grid_cache
#if (1)
int brgastro::grid_cache::_ra_grid_change_num_ = 0;
int brgastro::grid_cache::_dec_grid_change_num_ = 0;
int brgastro::grid_cache::_z_grid_change_num_ = 0;
BRG_ANGLE brgastro::grid_cache::_ra_grid_min_ = -pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_max_ = pi;
BRG_ANGLE brgastro::grid_cache::_ra_grid_step_ = pi / 8;
BRG_ANGLE brgastro::grid_cache::_dec_grid_min_ = -pi / 2;
BRG_ANGLE brgastro::grid_cache::dec_grid_max_val = pi / 2;
BRG_ANGLE brgastro::grid_cache::_dec_grid_step_ = pi / 8;
double brgastro::grid_cache::_z_grid_min_ = 0;
double brgastro::grid_cache::_z_grid_max_ = 2;
double brgastro::grid_cache::_z_grid_step_ = 0.1;
#endif // End initialisation for brgastro::grid_cache

// Initialisation for brgastro::dfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( dfa_cache, 0, 5, 0.01 );

// Initialisation for brgastro::add_cache
#if (1)
bool brgastro::add_cache::loaded = false;
std::string brgastro::add_cache::file_name = "brgastro_add_cache.dat";
double brgastro::add_cache::z_min = 0; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::add_cache::z_max = 5;
double brgastro::add_cache::z_step = 0.01;
int brgastro::add_cache::z_res = (int)( ( add_cache::z_max - add_cache::z_min )
		/ add_cache::z_step ) + 1;
std::string brgastro::add_cache::header_string = "# ADD cache v1.0";
int brgastro::add_cache::sig_digits = 8;
std::vector< std::vector< double > > brgastro::add_cache::dmod;
#endif // End initialization for brgastro::add_cache

// Initialisation for brgastro::tfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( tfa_cache, 0.001, 1.02, 0.001 );

// Initialisation for brgastro::tNFW_sig_cache
#if (1)
bool brgastro::tNFW_sig_cache::loaded = false;
std::string brgastro::tNFW_sig_cache::file_name = "brgastro_NFW_sig_cache.dat";
double brgastro::tNFW_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_sig_cache::halo_z_res = (int)( ( tNFW_sig_cache::halo_z_max
		- tNFW_sig_cache::halo_z_min ) / tNFW_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_sig_cache::halo_m_min = 10;
double brgastro::tNFW_sig_cache::halo_m_max = 16;
double brgastro::tNFW_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_sig_cache::halo_m_res = (int)( ( tNFW_sig_cache::halo_m_max
		- tNFW_sig_cache::halo_m_min ) / tNFW_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_sig_cache::r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
std::string brgastro::tNFW_sig_cache::header_string = "# NFW_sig v1.0";
int brgastro::tNFW_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< double > > > brgastro::tNFW_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_sig_cache

// Initialization for brgastro::tNFW_offset_sig_cache
#if (1)
bool brgastro::tNFW_offset_sig_cache::loaded = false;
std::string brgastro::tNFW_offset_sig_cache::file_name =
		"brgastro_NFW_offset_sig_cache.dat";
double brgastro::tNFW_offset_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_offset_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_offset_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_offset_sig_cache::halo_z_res =
		(int)( ( tNFW_offset_sig_cache::halo_z_max
				- tNFW_offset_sig_cache::halo_z_min )
				/ tNFW_offset_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_offset_sig_cache::halo_m_min = 10;
double brgastro::tNFW_offset_sig_cache::halo_m_max = 16;
double brgastro::tNFW_offset_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_offset_sig_cache::halo_m_res =
		(int)( ( tNFW_offset_sig_cache::halo_m_max
				- tNFW_offset_sig_cache::halo_m_min )
				/ tNFW_offset_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_offset_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_offset_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_offset_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_offset_sig_cache::r_res =
		(int)( ( tNFW_offset_sig_cache::r_max - tNFW_offset_sig_cache::r_min )
				/ tNFW_offset_sig_cache::r_step ) + 1;
double brgastro::tNFW_offset_sig_cache::offset_r_min = -1;
double brgastro::tNFW_offset_sig_cache::offset_r_max = 4;
double brgastro::tNFW_offset_sig_cache::offset_r_step = 0.1;
int brgastro::tNFW_offset_sig_cache::offset_r_res =
		(int)( ( tNFW_offset_sig_cache::offset_r_max
				- tNFW_offset_sig_cache::offset_r_min )
				/ tNFW_offset_sig_cache::offset_r_step ) + 1;
std::string brgastro::tNFW_offset_sig_cache::header_string =
		"# NFW_offset_sig v1.0";
int brgastro::tNFW_offset_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< std::vector< double > > > > brgastro::tNFW_offset_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_offset_sig_cache

// Initialization for brgastro::tNFW_group_sig_cache
#if (1)
bool brgastro::tNFW_group_sig_cache::loaded = false;
std::string brgastro::tNFW_group_sig_cache::file_name =
		"brgastro::tNFW_group_sig_cache.dat";
double brgastro::tNFW_group_sig_cache::halo_z_min = 0.1; // These values will only be used if the cache needs to be created and they aren't overridden
double brgastro::tNFW_group_sig_cache::halo_z_max = 1.3;
double brgastro::tNFW_group_sig_cache::halo_z_step = 0.1;
int brgastro::tNFW_group_sig_cache::halo_z_res =
		(int)( ( tNFW_group_sig_cache::halo_z_max
				- tNFW_group_sig_cache::halo_z_min )
				/ tNFW_group_sig_cache::halo_z_step ) + 1;
double brgastro::tNFW_group_sig_cache::halo_m_min = 10;
double brgastro::tNFW_group_sig_cache::halo_m_max = 16;
double brgastro::tNFW_group_sig_cache::halo_m_step = 0.01;
int brgastro::tNFW_group_sig_cache::halo_m_res =
		(int)( ( tNFW_group_sig_cache::halo_m_max
				- tNFW_group_sig_cache::halo_m_min )
				/ tNFW_group_sig_cache::halo_m_step ) + 1;
double brgastro::tNFW_group_sig_cache::r_min = 1 * unitconv::kpctom;
double brgastro::tNFW_group_sig_cache::r_max = 5000 * unitconv::kpctom;
double brgastro::tNFW_group_sig_cache::r_step = 1 * unitconv::kpctom;
int brgastro::tNFW_group_sig_cache::r_res =
		(int)( ( tNFW_group_sig_cache::r_max - tNFW_group_sig_cache::r_min )
				/ tNFW_group_sig_cache::r_step ) + 1;
double brgastro::tNFW_group_sig_cache::group_c_min = 2;
double brgastro::tNFW_group_sig_cache::group_c_max = 8;
double brgastro::tNFW_group_sig_cache::group_c_step = 1;
int brgastro::tNFW_group_sig_cache::group_c_res =
		(int)( ( tNFW_group_sig_cache::group_c_max
				- tNFW_group_sig_cache::group_c_min )
				/ tNFW_group_sig_cache::group_c_step ) + 1;
std::string brgastro::tNFW_group_sig_cache::header_string =
		"# NFW_group_sig v1.0";
int brgastro::tNFW_group_sig_cache::sig_digits = 8;
std::vector< std::vector< std::vector< std::vector< double > > > > brgastro::tNFW_group_sig_cache::signal;
#endif // End initialization for brgastro::tNFW_group_sig_cache

#endif // end Static Class Initialisation

/** Class Method Definitions **/
#if (1)
// brgastro::redshift_obj class methods
#if (1)
const int brgastro::redshift_obj::z_grid() const
{
	if ( _z_grid_cached_ )
	{
		if ( _local_z_grid_change_num_ == grid_cache().z_grid_change_num() )
			return _z_grid_;
	}
	_z_grid_ = brgastro::get_z_grid( z() );
	_z_grid_cached_ = true;
	_local_z_grid_change_num_ = grid_cache().z_grid_change_num();
	return _z_grid_;
}

const BRG_UNITS brgastro::redshift_obj::H( const double init_test_z ) const
{
	double test_z = init_test_z;
	if ( test_z == -1 ) // This means we're taking the default argument here, so use redshift of the object
	{
		test_z = z();
	}
	// Friedmann equation, assuming omega = -1
	return H_0
			* sqrt(
					Omega_r * std::pow( 1 + test_z, 4 )
							+ Omega_m * std::pow( 1 + test_z, 3 )
							+ Omega_k * pow( 1 + test_z, 2 ) + Omega_l );
}

#endif

// brgastro::sky_obj class methods
#if (1)

brgastro::sky_obj::sky_obj( const BRG_ANGLE init_ra, const BRG_ANGLE init_dec,
		const double init_z, const BRG_ANGLE init_ra_err,
		const BRG_ANGLE init_dec_err, const double init_z_err )
{
	partial_clear();
	set_ra_dec_z_err( init_ra, init_dec, init_z, init_ra_err, init_dec_err,
			init_z_err );
}

const int brgastro::sky_obj::clear()
{
	set_ra_dec_z_err( 0, 0, 0, 0, 0, 0 );
	return partial_clear();
}

const int brgastro::sky_obj::partial_clear()
{
	_index_ = 0;
	_ID_ = "0";
	_weight_ = 1;
	_ra_grid_ = _dec_grid_ = 0;
	_ra_grid_cached_ = _dec_grid_cached_ = false;
	_local_ra_grid_change_num_ = -1;
	_local_dec_grid_change_num_ = -1;
	return 0;
}

const int brgastro::sky_obj::set_ra( const BRG_ANGLE &new_ra )
{
	if ( _ra_ == new_ra )
		return 0;
	_ra_ = new_ra;
	_ra_grid_cached_ = false;
	return 0;
}
const int brgastro::sky_obj::set_ra_err( const BRG_ANGLE &new_ra_err )
{
	_ra_err_ = new_ra_err;
	return 0;
}
const int brgastro::sky_obj::set_dec( const BRG_ANGLE &new_dec )
{
	if ( _dec_ == new_dec )
		return 0;
	_dec_ = new_dec;
	_dec_grid_cached_ = false;
	return 0;
}
const int brgastro::sky_obj::set_dec_err( const BRG_ANGLE &new_dec_err )
{
	_dec_err_ = new_dec_err;
	return 0;
}
const int brgastro::sky_obj::set_ra_dec( const BRG_ANGLE &new_ra,
		const BRG_ANGLE &new_dec )
{
	if ( set_ra( new_ra ) )
		return errorNOS();
	return set_dec( new_dec );
}
const int brgastro::sky_obj::set_ra_dec_z( const BRG_ANGLE &new_ra,
		const BRG_ANGLE &new_dec, const double new_z )
{
	if ( set_ra( new_ra ) )
		return errorNOS();
	if ( set_dec( new_dec ) )
		return errorNOS();
	return set_z( new_z );
}
const int brgastro::sky_obj::set_ra_dec_z_err( const BRG_ANGLE &new_ra,
		const BRG_ANGLE &new_dec, const double new_z,
		const BRG_ANGLE &new_ra_err, const BRG_ANGLE &new_dec_err,
		const double new_z_err )
{
	if ( set_ra( new_ra ) )
		return errorNOS();
	if ( set_dec( new_dec ) )
		return errorNOS();
	if ( set_z( new_z ) )
		return errorNOS();
	if ( set_ra_err( new_ra_err ) )
		return errorNOS();
	if ( set_dec_err( new_dec_err ) )
		return errorNOS();
	return set_z_err( new_z_err );
}
const int brgastro::sky_obj::set_ra_dec_err( const BRG_ANGLE &new_ra,
		const BRG_ANGLE &new_dec, const BRG_ANGLE &new_ra_err,
		const BRG_ANGLE &new_dec_err )
{
	if ( set_ra( new_ra ) )
		return errorNOS();
	if ( set_dec( new_dec ) )
		return errorNOS();
	if ( set_ra_err( new_ra_err ) )
		return errorNOS();
	return set_dec_err( new_dec_err );
}
const int brgastro::sky_obj::set_weight( const double new_weight )
{
	_weight_ = new_weight;
	return 0;
}
const int brgastro::sky_obj::set_index( const int new_index )
{
	_index_ = new_index;
	return 0;
}
const int brgastro::sky_obj::set_ID( const std::string &new_ID )
{
	_ID_ = new_ID;
	return 0;
}

const int brgastro::sky_obj::ra_grid() const
{
	if ( _ra_grid_cached_ )
	{
		if ( _local_ra_grid_change_num_ == grid_cache().ra_grid_change_num() )
			return _ra_grid_;
	}
	_ra_grid_ = brgastro::get_ra_grid( ra() );
	_ra_grid_cached_ = true;
	_local_ra_grid_change_num_ = grid_cache().ra_grid_change_num();
	return _ra_grid_;
}
const int brgastro::sky_obj::dec_grid() const
{
	if ( _dec_grid_cached_ )
	{
		if ( _local_dec_grid_change_num_
				== grid_cache().dec_grid_change_num() )
			return _dec_grid_;
	}
	_dec_grid_ = brgastro::get_dec_grid( dec() );
	_dec_grid_cached_ = true;
	_local_dec_grid_change_num_ = grid_cache().dec_grid_change_num();
	return _dec_grid_;
}

#endif // end brgastro::sky_obj class functions

// brgastro::galaxy class methods
#if (1)
brgastro::galaxy::galaxy()
{

	clear();

}

const int brgastro::galaxy::clear()
{

	stellar_mass = 0;
	umag = gmag = rmag = imag = zmag = 0;
	umag_err = gmag_err = rmag_err = imag_err = zmag_err = 0;

	z_phot = z();
	z_phot_err = 0;
	odds = 1;
	phot_template = 0;

	host_group = 0;
	host_group_index = -1;

	return ( sky_obj::clear() );
}

#endif // end brgastro::galaxy class functions

// brgastro::group class methods
#if (1)

// Constructor
brgastro::group::group( int init_num_members )
{
	if ( init_num_members < 0 )
		init_num_members = 0;

	z_phot = z_phot_err = 0;
	odds = 1;

	BCG_index = 0;

	num_members = init_num_members;
	if ( num_members > 0 )
	{
		member_indices.reserve( num_members );
		members.reserve( num_members );
	}
	for ( int i = 0; i < num_members; i++ )
	{
		member_indices[i] = -1;
		members[i] = 0;
	}
}

// Copy constructor
brgastro::group::group( const group &other_group ) :
		sky_obj( other_group )
{

	z_phot = other_group.z_phot;
	z_phot_err = other_group.z_phot_err;

	odds = other_group.odds;

	BCG_index = other_group.BCG_index;

	num_members = other_group.num_members;
	if ( num_members > 0 )
	{
		member_indices.reserve( num_members );
		members.reserve( num_members );
	}
	for ( int i = 0; i < num_members; i++ )
	{
		member_indices[i] = other_group.member_indices[i];
		members[i] = other_group.members[i];
	}

}

// Destructor
brgastro::group::~group()
{
}

// Other functions
const int brgastro::group::clear()
{
	group();
	return 0;
}

const int brgastro::group::resize( int new_num_members )
{
	member_indices.resize( new_num_members, -1 );
	members.resize( new_num_members, 0 );
	num_members = new_num_members;
	return 0;
}

const int brgastro::group::set_member( int member_index, galaxy * new_member,
		const bool silent )
{
	if ( member_index >= num_members )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Member index out of bounds in group::set_member\n";
		return OUT_OF_BOUNDS_ERROR;
	}
	members[member_index] = new_member;
	member_indices[member_index] = new_member->index();
	return 0;
}
const int brgastro::group::set_member_index( int member_index,
		int new_member_index, const bool silent )
{
	if ( member_index >= num_members )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Member index out of bounds in group::set_member_index\n";
		return OUT_OF_BOUNDS_ERROR;
	}
	member_indices[member_index] = new_member_index;
	return 0;
}
const int brgastro::group::add_member( galaxy * new_member, const bool silent )
{
	int new_num_members = num_members + 1;
	resize( new_num_members );
	members[new_num_members - 1] = new_member;
	member_indices[new_num_members - 1] = new_member->index();
	return 0;
}
const int brgastro::group::remove_member( galaxy * rem_member,
		const bool silent )
{
	int i;
	for ( i = 0; i < num_members; i++ )
	{
		if ( members[i]->ID() == rem_member->ID() )
			break;
	}
	if ( i == num_members ) // Then we didn't find it
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not find member to remove from group.\n";
		return OUT_OF_BOUNDS_ERROR;
	}
	for ( int j = i; j < num_members - 1; j++ )
	{
		members[j] = members[j + 1];
		member_indices[j] = member_indices[j + 1];
	}
	resize( num_members - 1 );
	return 0;
}

#endif // End brgastro::group class functions

// brgastro::density_profile class methods
#if (1)
const BRG_UNITS brgastro::density_profile::Daccel( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_DISTANCE dr;
	BRG_UNITS a1, a2;
	// It's simplest, and a ton faster to just manually differentiate here.
	dr = max( r * SMALL_FACTOR, SMALL_FACTOR );
	a1 = accel( r, silent );
	a2 = accel( r + dr, silent );
	return ( a2 - a1 ) / safe_d( dr );
}

const BRG_DISTANCE brgastro::density_profile::rhmtot( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmtot_cached )
		return _rhmtot_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmtot();
	solve_rhm_function func( this, target_mass );
	solve_rhm_function *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmtot_cached = true;
		return _rhmtot_cache_ = 0;
	}

	if ( solve_grid( func_ptr, 1, 0., max_r, 10, 0., rhm_test ) ) // If we can't solve it
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius. Assuming it's zero.\n";

		_rhmtot_cache_ = 0;
		hmtot_cached = true;
		return _rhmtot_cache_;
	}
	else
	{
		_rhmtot_cache_ = std::fabs( rhm_test );
		hmtot_cached = true;
	}
	return _rhmtot_cache_;
}

const BRG_DISTANCE brgastro::density_profile::rhmvir( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmvir_cached )
		return _rhmvir_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmvir();
	solve_rhm_function func( this, target_mass );
	solve_rhm_function *func_ptr = &func;

	BRG_UNITS max_r( default_tau_factor * rvir() );
	BRG_UNITS rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmvir_cached = true;
		return _rhmvir_cache_ = 0;
	}

	if ( solve_grid( func_ptr, 1, 0., max_r, 10, 0., rhm_test ) ) // If we can't solve it
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius.\n";
		return -1;
	}
	else
	{
		_rhmvir_cache_ = max(0,std::fabs( rhm_test ));
		hmvir_cached = true;
	}
	return _rhmvir_cache_;
}

const BRG_MASS brgastro::density_profile::enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	if ( r == 0 )
		return 0;
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::spherical_density_function func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( r_to_use ), out_params( 0 );
	if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate enclosed mass of density profile.\n";
	}
	return out_params;
}

const BRG_UNITS brgastro::density_profile::proj_dens( const BRG_DISTANCE &R,
		const bool silent ) const
{
	BRG_DISTANCE R_to_use = std::fabs( R );
	double inf_factor = 20;
	brgastro::cylindrical_density_function func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( inf_factor * rvir() ),
			out_params( 0 );
	if ( R_to_use == 0 )
	{
		// In this case, we might be integrating over a singularity, so the trapezoid method is safer
		const int num_steps = 10000;

		if ( brgastro::integrate( &func, num_in_params, min_in_params,
				max_in_params, num_steps, num_out_params, out_params ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not integrate projected density of density profile.\n";
		}
	}
	else
	{
		if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
				max_in_params, num_out_params, out_params ) )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not integrate projected density of density profile.\n";
		}
	}
	return 2 * out_params;
}

const BRG_MASS brgastro::density_profile::proj_enc_mass( const BRG_DISTANCE &R,
		const bool silent ) const
{
	if ( R == 0 )
		return 0;
	BRG_DISTANCE R_to_use = std::fabs( R );
	brgastro::cylindrical_density_function func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( 0 ), max_in_params( R_to_use ), out_params( 0 );
	if ( brgastro::integrate_Rhomberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate projected enclosed mass of density profile.\n";
	}
	return out_params;
}

const BRG_UNITS brgastro::density_profile::offset_WLsig( const BRG_DISTANCE &R,
		const BRG_DISTANCE &offset_R, const bool silent ) const
{
	unsigned int num_out_params = 1;
	if ( offset_R == 0 )
		return deltasigma( R );
	BRG_DISTANCE R_to_use = std::fabs( R );
	BRG_DISTANCE offset_R_to_use = std::fabs( R );
	offset_ring_dens_function ringfunc( this, offset_R_to_use, R_to_use );
	offset_circ_dens_function circfunc( this, offset_R_to_use );

	BRG_UNITS min_in_params_ring( 0 );
	BRG_UNITS max_in_params_ring( 0 );
	BRG_UNITS out_params_ring( 0 );
	std::vector< BRG_UNITS > min_in_params_circ( 2, 0 );
	std::vector< BRG_UNITS > max_in_params_circ( 2, 0 );
	std::vector< BRG_UNITS > out_params_circ( 1, 0 );
	BRG_UNITS circmean;
	BRG_UNITS ringmean;
	BRG_UNITS result;

	double precision = 0.001;

	min_in_params_ring = 0;
	max_in_params_ring = pi; // We'll double it later, due to symmetry

	min_in_params_circ[0] = R_to_use * precision;
	max_in_params_circ[0] = R_to_use;
	min_in_params_circ[1] = 0;
	max_in_params_circ[1] = pi; // We'll double it later, due to symmetry

	if ( brgastro::integrate_Rhomberg( &ringfunc, 1, min_in_params_ring,
			max_in_params_ring, num_out_params, out_params_ring, precision ) )
	{
		if ( !silent )
			std::cerr << "ERROR: Could not integrate in offset_WLsig" << endl;
		return 0;
	}

	ringmean = 2 * out_params_ring / pi;

	if ( brgastro::integrate_Rhomberg( &circfunc, 2, min_in_params_circ,
			max_in_params_circ, num_out_params, out_params_circ, precision ) )

	{
		if ( !silent )
			std::cerr
					<< "ERROR: Could not integrate in unitless_NFW_offset_WLsig"
					<< endl;
		return 0;
	}

	circmean = out_params_circ[0] / ( pi * R_to_use );

	result = circmean - ringmean;

	return result;
}
const BRG_UNITS brgastro::density_profile::group_WLsig( const BRG_DISTANCE &r,
		const double group_c, const bool silent ) const
{
	BRG_DISTANCE r_to_use = std::fabs( r );
	brgastro::offset_WLsig_function func( this, r_to_use );
	brgastro::offset_WLsig_weight_function weight_func( this, group_c );
	unsigned int num_in_params = 1, num_out_params = 0;
	BRG_UNITS min_in_params( SMALL_FACTOR ), max_in_params( rvir() ),
			out_params( 0 );
	if ( brgastro::integrate_weighted_Rhomberg( &func, &weight_func,
			num_in_params, min_in_params, max_in_params, num_out_params,
			out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate group deltasigma of density profile.\n";
	}
	return out_params;
}

#endif

// brgastro::tNFW_profile class methods
#if (1)

#if (1) // Constructors
brgastro::tNFW_profile::tNFW_profile()
{
	_mvir0_ = 0;
	_c_ = 0;
	_tau_ = 0;
}

brgastro::tNFW_profile::tNFW_profile( const BRG_MASS &init_mvir0,
		const double init_z, const double init_c, const double init_tau ) :
		brgastro::redshift_obj( init_z )
{
	_mvir0_ = init_mvir0;
	if ( init_c <= 0 )
	{
		_c_ = cfm( init_mvir0, init_z );
	}
	else
	{
		_c_ = init_c;
	}
	if ( init_tau < 0 )
	{
		_tau_ = default_tau_factor * _c_;
	}
	else
	{
		_tau_ = init_tau;
	}
}

#endif // End constructors

// Destructor
brgastro::tNFW_profile::~tNFW_profile()
{
}

#if (1) // Set functions

const int brgastro::tNFW_profile::set_mvir( const BRG_MASS &new_halo_mass,
		const bool silent )
{
	_mvir0_ = new_halo_mass;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_tau( const double new_halo_tau,
		const bool silent )
{
	_tau_ = new_halo_tau;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_c( const double new_halo_c,
		const bool silent )
{
	_c_ = new_halo_c;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_z( const double new_z )
{
	redshift_obj::set_z( new_z );
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
const int brgastro::tNFW_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< BRG_UNITS > &parameters, const bool silent )
{
	if ( ( num_parameters != 4 ) || ( num_parameters != parameters.size() ) )
	{
		if ( !silent )
			if ( !silent )
				std::cerr
						<< "ERROR: Invalid number of parameters passed to tNFW_profile::set_parameters.\n"
						<< "Four are required for both num_parameters and parameters.size().\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( set_mvir( parameters.at( 0 ) ) )
		return errorNOS( silent );
	if ( set_z( parameters.at( 1 ) ) )
		return errorNOS( silent );
	if ( parameters.at( 2 ) <= 0 )
	{
		if ( set_c( cfm( parameters.at( 0 ), parameters.at( 1 ) ) ) )
			return errorNOS( silent );
	}
	else
	{
		if ( set_c( parameters.at( 2 ) ) )
			return errorNOS( silent );
	}
	if ( parameters.at( 3 ) <= 0 )
	{
		if ( set_tau( default_tau_factor * _c_ ) )
			return errorNOS( silent );
	}
	else
	{
		if ( set_tau( parameters.at( 3 ) ) )
			return errorNOS( silent );
	}
	return 0;
}

#endif // end set functions

const BRG_MASS brgastro::tNFW_profile::mvir() const
{
	return _mvir0_; // Not technically correct, but close enough for our purposes
}
const BRG_MASS brgastro::tNFW_profile::mvir0() const
{
	return _mvir0_;
}

const double brgastro::tNFW_profile::tau() const
{
	return _tau_;
}
const double brgastro::tNFW_profile::c() const
{
	return _c_;
}

const BRG_MASS brgastro::tNFW_profile::mtot() const
{
	return _mvir0_ * brgastro::mftau( _tau_, _c_ );
}

const BRG_VELOCITY brgastro::tNFW_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
const BRG_DISTANCE brgastro::tNFW_profile::rvir() const
{
	return vvir() / H() / 10;
}
const BRG_DISTANCE brgastro::tNFW_profile::rs() const
{
	return rvir() / _c_;
}
const BRG_DISTANCE brgastro::tNFW_profile::rt( const bool silent ) const
{
	return rvir() / _tau_;
}

const BRG_MASS brgastro::tNFW_profile::hmvir() const
{
	return enc_mass( rvir() ) / 2;
}

const BRG_UNITS brgastro::tNFW_profile::quick_WLsig( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result;

	result = brgastro::tNFW_sig_cache().get( z(), _mvir0_, r, silent );
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::quick_offset_WLsig(
		const BRG_DISTANCE &r, const BRG_DISTANCE &offset_r,
		const bool silent ) const
{
	BRG_UNITS result;
	result = brgastro::tNFW_offset_sig_cache().get( z(), _mvir0_, r, offset_r,
			silent );
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::semiquick_group_WLsig(
		const BRG_DISTANCE &r, const double group_c, const bool silent ) const
{
	BRG_UNITS result;
	if ( !silent )
		std::cerr << "ERROR: Placeholder function being used.\n";
	result = 0; // Placeholder!!!
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::quick_group_WLsig(
		const BRG_DISTANCE &r, const double group_c, const bool silent ) const
{
	BRG_UNITS result;
	result = brgastro::tNFW_group_sig_cache().get( z(), _mvir0_, r, group_c,
			silent );
	return result;
}

const BRG_UNITS brgastro::tNFW_profile::dens( const BRG_DISTANCE &r ) const
{
	BRG_UNITS result, rho_c;

	double d_c, x, tau_use;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;
	d_c = delta_c( _c_ );
	rho_c = 3 * H() * H() / ( 8 * pi * Gc );
	x = r / rs();

	result = ( d_c * rho_c ) / ( x * std::pow( 1 + x, 2 ) )
			* std::pow( tau_use, 2 )
			/ ( std::pow( tau_use, 2 ) + std::pow( x, 2 ) );

	return result;
}
const BRG_UNITS brgastro::tNFW_profile::proj_dens( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c;
	double d_c, x, tau_use, fx, lx;

	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;

	d_c = delta_c( _c_ );
	rho_c = 3 * H() * H() / ( 8 * pi * Gc );
	x = r / rs();

	if ( x == 1 )
		fx = 1;
	else if ( x > 1 )
		fx = acos( 1 / x ) / sqrt( x * x - 1 );
	else
		fx = -log( 1 / x - sqrt( 1 / ( x * x ) - 1 ) ) / sqrt( 1 - x * x );
	lx = log( x / ( sqrt( tau_use * tau_use + x * x ) + tau_use ) );
	if ( x == 1 )
		result =
				( 4 * pi * rs() * d_c * rho_c ) * tau_use * tau_use
						/ ( 2 * pi * std::pow( tau_use * tau_use + 1, 2 ) )
						* ( ( tau_use * tau_use + 1 ) / 3 + 2 * fx
								- pi / sqrt( tau_use * tau_use + x * x )
								+ ( tau_use * tau_use - 1 ) * lx
										/ ( tau_use
												* sqrt(
														tau_use * tau_use
																+ x * x ) ) );
	else
		result =
				( 4 * pi * rs() * d_c * rho_c ) * tau_use * tau_use
						/ ( 2 * pi * std::pow( tau_use * tau_use + 1, 2 ) )
						* ( ( tau_use * tau_use + 1 ) / ( x * x - 1 )
								* ( 1 - fx ) + 2 * fx
								- pi / sqrt( tau_use * tau_use + x * x )
								+ ( tau_use * tau_use - 1 ) * lx
										/ ( tau_use
												* sqrt(
														tau_use * tau_use
																+ x * x ) ) );
	return result;
}
const BRG_MASS brgastro::tNFW_profile::enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c;
	BRG_MASS m0;
	// Result here integrated with Wolfram Alpha
	double d_c, x, tau_use;
	if ( _tau_ < 0 )
		tau_use = default_tau_factor * _c_;
	else if (_tau_ == 0)
		return 0;
	else
		tau_use = _tau_;

	d_c = delta_c( _c_ );
	rho_c = 3 * H() * H() / ( 8 * pi * Gc );
	x = r / rs();

	m0 = ( std::pow( rs(), 3 ) * d_c * rho_c )
			* ( 2 * pi * std::pow( tau_use, 2 )
					* ( 2 * ( 1 + std::pow( tau_use, 2 ) )
							- ( -1 + std::pow( tau_use, 2 ) ) * 2
									* log( tau_use ) ) )
			/ std::pow( 1 + std::pow( tau_use, 2 ), 2 );
	result = ( std::pow( rs(), 3 ) * d_c * rho_c )
			* ( 2 * pi * std::pow( tau_use, 2 )
					* ( ( 2 * ( 1 + std::pow( tau_use, 2 ) ) ) / ( 1 + x )
							+ 4 * tau_use * atan( x / tau_use )
							+ 2 * ( -1 + std::pow( tau_use, 2 ) )
									* log( 1 + x )
							- ( -1 + std::pow( tau_use, 2 ) )
									* log(
											std::pow( tau_use, 2 )
													+ std::pow( x, 2 ) ) ) )
			/ std::pow( 1 + std::pow( tau_use, 2 ), 2 ) - m0;
	return result;
}
const BRG_UNITS brgastro::tNFW_profile::proj_enc_dens( const BRG_DISTANCE &r,
		const bool silent ) const
{
	BRG_UNITS result, rho_c;
	//Takes M in kg, r in kpc
	double d_c, x, tau_use, fx, lx;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;

	d_c = delta_c( _c_ );
	rho_c = 3 * H() * H() / ( 8 * pi * Gc );
	x = r / rs();
	if ( x == 1 )
		fx = 1;
	else if ( x > 1 )
		fx = acos( 1 / x ) / sqrt( x * x - 1 );
	else
		fx = -log( 1 / x - sqrt( 1 / ( x * x ) - 1 ) ) / sqrt( 1 - x * x );
	lx = log( x / ( sqrt( tau_use * tau_use + x * x ) + tau_use ) );

	result =
			( 4 * pi * rs() * d_c * rho_c ) * tau_use * tau_use
					/ ( pi * x * x * std::pow( tau_use * tau_use + 1, 2 ) )
					* ( ( tau_use * tau_use + 1 + 2 * ( x * x - 1 ) ) * fx
							+ pi * tau_use
							+ ( tau_use * tau_use - 1 ) * log( tau_use )
							+ sqrt( tau_use * tau_use + x * x )
									* ( lx * ( tau_use * tau_use - 1 )
											/ tau_use - pi ) );
	return result;
}
const BRG_MASS brgastro::tNFW_profile::proj_enc_mass( const BRG_DISTANCE &r,
		const bool silent ) const
{
	return proj_enc_dens( r, silent ) * pi * r * r;
}

const int brgastro::tNFW_profile::get_parameters( int & num_parameters,
		std::vector< BRG_UNITS > & parameters, const bool silent ) const
{
	num_parameters = 4;
	parameters.resize( num_parameters );

	try
	{
		parameters.at( 0 ) = _mvir0_;
		parameters.at( 1 ) = z();
		parameters.at( 2 ) = _c_;
		parameters.at( 3 ) = _tau_;
	}
	catch ( std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}

const int brgastro::tNFW_profile::get_parameter_names( int & num_parameters,
		std::vector< std::string > & parameter_names, const bool silent ) const
{
	num_parameters = 4;
	parameter_names.resize( num_parameters );

	try
	{
		parameter_names.at( 0 ) = "mvir0";
		parameter_names.at( 1 ) = "z";
		parameter_names.at( 2 ) = "c";
		parameter_names.at( 3 ) = "tau";
	}
	catch ( std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}

const int brgastro::tNFW_profile::truncate_to_fraction( const double f,
		const bool silent )
{
	if ( f <= 0 )
	{
		if ( f < 0 )
			if ( !silent )
				std::cerr
						<< "WARNING: Cannot truncate to negative fraction. Truncating to zero instead.\n";
		_tau_ = 0;
		override_rhmvir( 0 );
		override_rhmtot( 0 );
	}
	else
	{
		double new_tau_val = brgastro::taufm( f, _c_, _tau_ );
		if ( new_tau_val < 0 )
		{
			if ( !silent )
				std::cerr
						<< "WARNING: Could not solve for new tau value in truncate_to_fraction.\n"
						<< "tau will remain unchanged.\n";
		}
		else
		{
			_tau_ = new_tau_val;
		}
		hmvir_cached = false;
		hmtot_cached = false;
	}
	return 0;

}

#endif // end tNFW_profile functions

// brgastro::point_mass_profile class methods
#if (1)

#if (1) // Constructors
brgastro::point_mass_profile::point_mass_profile()
{
	_mass_ = 0;
}

brgastro::point_mass_profile::point_mass_profile( const BRG_MASS init_mass,
		const double init_z )
{
	set_mvir( init_mass );
	set_z( init_z );
}

#endif // End constructors

// Destructor
brgastro::point_mass_profile::~point_mass_profile()
{
}

#if (1) // Set functions

const int brgastro::point_mass_profile::set_mvir(
		const BRG_MASS &new_halo_mass, bool silent )
{
	_mass_ = new_halo_mass;
	return 0;
}
const int brgastro::point_mass_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< BRG_UNITS > &parameters, bool silent )
{

	if ( ( num_parameters != 2 ) || ( num_parameters != parameters.size() ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Invalid number of parameters passed to tNFW_profile::set_parameters.\n"
					<< "Four are required for both num_parameters and parameters.size().\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	if ( set_mvir( parameters.at( 0 ) ) )
		return errorNOS();
	if ( set_z( parameters.at( 1 ) ) )
		return errorNOS();
	return 0;
}

#endif // end set functions

const BRG_MASS brgastro::point_mass_profile::mvir() const
{
	return _mass_;
}
const BRG_MASS brgastro::point_mass_profile::mass() const
{
	return _mass_;
}

const BRG_MASS brgastro::point_mass_profile::mtot() const
{
	return _mass_;
}

const BRG_VELOCITY brgastro::point_mass_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
const BRG_DISTANCE brgastro::point_mass_profile::rvir() const
{
	return vvir() / H() / 10;
}
const BRG_DISTANCE brgastro::point_mass_profile::rs() const
{
	return 0;
}
const BRG_DISTANCE brgastro::point_mass_profile::rt(
		const bool silent) const
{
	return 0;
}

const BRG_MASS brgastro::point_mass_profile::hmvir() const
{
	return 0;
}

const BRG_UNITS brgastro::point_mass_profile::dens(
		const BRG_DISTANCE &r ) const
{
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-3,0,1,0,0,0);
#else
	double result = 0;
#endif

	result = ( r == 0 ? DBL_MAX : 0 );

	return result;
}
const BRG_UNITS brgastro::point_mass_profile::proj_dens(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-2,0,1,0,0,0);
#else
	double result = 0;
#endif

	result = ( r == 0 ? DBL_MAX : 0 );

	return result;
}
const BRG_UNITS brgastro::point_mass_profile::enc_dens(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return enc_mass( r ) / ( 4. / 3. * pi * std::pow( r, 3 ) );
}
const BRG_MASS brgastro::point_mass_profile::enc_mass(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return _mass_;
}
const BRG_UNITS brgastro::point_mass_profile::proj_enc_dens(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return proj_enc_mass( r ) / ( pi * r * r );
}
const BRG_MASS brgastro::point_mass_profile::proj_enc_mass(
		const BRG_DISTANCE &r,
		const bool silent ) const
{
	return _mass_;
}

const int brgastro::point_mass_profile::get_parameters( int & num_parameters,
		std::vector< BRG_UNITS > & parameters, const bool silent ) const
{
	num_parameters = 2;
	parameters.resize( num_parameters );

	try
	{
		parameters.at( 0 ) = _mass_;
		parameters.at( 1 ) = z();
	}
	catch ( std::exception & )
	{
		return errorNOS();
	}
	return 0;
}

const int brgastro::point_mass_profile::get_parameter_names(
		int & num_parameters, std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	num_parameters = 2;
	parameter_names.resize( num_parameters );

	try
	{
		parameter_names.at( 0 ) = "mass";
		parameter_names.at( 1 ) = "z";
	}
	catch ( std::exception & )
	{
		return errorNOS();
	}
	return 0;
}

const int brgastro::point_mass_profile::truncate_to_fraction( const double f,
		const bool silent )
{
	if ( ( f <= 0 ) || ( isbad( f ) ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Bad, negative, or zero f passed to truncate_to_fraction.\n"
					<< "Will truncate to zero instead.\n";
		_mass_ = 0;
	}
	else
	{
		_mass_ *= f;
	}
	return 0;

}

#endif // end point_mass profile functions

// brgastro::dfa_cache class methods
#if (1)
const int brgastro::dfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = brgastro::integrate_add( 0, in_params );
	}
	catch(std::exception &e)
	{
		return LOWER_LEVEL_ERROR;
	}
	catch(...)
	{
		return LOWER_LEVEL_ERROR;
	}

	return 0;
}

#endif // end brgastro::dfa_cache functions

// brgastro::add_cache class methods
#if (1)

const int brgastro::add_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j;

	if ( loaded )
		return 0;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::add_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> z_min >> z_max >> z_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		z_res = (int)( ( z_max - z_min ) / z_step ) + 1;
		make_array2d( dmod, z_res, z_res );

		// Read in data
		i = j = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			dmod[i][j] = temp_data;
			j++;
			if ( j >= z_res )
			{
				j = 0;
				i++;
				if ( i >= z_res )
				{
					break;
				}
			}
		}

		// Check that it was all read properly
		if ( i < z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::add_cache::unload()
{
	loaded = false;
	return del_array2d( dmod, z_res );
}

const int brgastro::add_cache::calc( const bool silent )
{
	int i = 0, j = 0;

	// Test that range is sane
	if ( ( z_max <= z_min ) || ( z_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to add_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	z_res = (int)( ( z_max - z_min ) / z_step ) + 1;
	if ( make_array2d( dmod, z_res, z_res ) )
		return 1;

	// Calculate data
	for ( double z1 = z_min; z1 < z_max; z1 += z_step )
	{
		for ( double z2 = z_min; z2 < z_max; z2 += z_step )
		{
			dmod[i][j] = brgastro::integrate_add( z1, z2 );
			j++;
		}
		j = 0;
		i++;
	}

	return 0;
}

const int brgastro::add_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name, true ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << z_min << "\t" << z_max << "\t" << z_step << "\n";

	// Output data
	for ( int i = 0; i < z_res; i++ )
	{
		for ( int j = 0; j < z_res; j++ )
		{
			if ( !( out_file << dmod[i][j] << "\t" ) )
				return errorNOS();
		}
		out_file << "\n";
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::add_cache::set_file_name( const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		return unload();
	}
	return 0;
}

const int brgastro::add_cache::set_range( const double new_z_min,
		const double new_z_max, const double new_z_step, const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( z_min != new_z_min ) || ( z_max != new_z_max )
			|| ( z_step != new_z_step ) )
	{
		z_min = new_z_min;
		z_max = new_z_max;
		z_step = new_z_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::add_cache::set_precision( const int new_precision,
		const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( max( new_precision, 0 ), DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr << "ERROR: Precision for add_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::add_cache::get( const double z1, const double z2,
		const bool silent )
{
	double z1lo, z1hi, z2lo, z2hi, z1wlo, z1whi, z2wlo, z2whi;
	int z1_i, z2_i; // Lower nearby array point
	BRG_DISTANCE result = 0;

	if ( !loaded )
	{
		#pragma omp critical(add_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Could neither load nor calculate add_cache!\n";
			return result;
		}
	}

	z1_i = (int)( ( z1 - z_min ) / z_step );
	z1_i = max( z1_i, 0 );
	z1_i = min( z1_i, z_res - 2 );

	z1lo = z_min + z_step * z1_i;
	z1hi = z_min + z_step * ( z1_i + 1 );
	z1wlo = z1 - z1lo;
	z1whi = z1hi - z1;

	z2_i = (int)( ( z2 - z_min ) / z_step );
	z2_i = max( z2_i, 0 );
	z2_i = min( z2_i, z_res - 2 );

	z2lo = z_min + z_step * z2_i;
	z2hi = z_min + z_step * ( z2_i + 1 );
	z2wlo = z2 - z2lo;
	z2whi = z2hi - z2;

	double totweight = z_step * z_step;

	result += dmod[z1_i][z2_i] * z1whi * z2whi;
	result += dmod[z1_i][z2_i + 1] * z1whi * z2wlo;
	result += dmod[z1_i + 1][z2_i] * z1wlo * z2whi;
	result += dmod[z1_i + 1][z2_i + 1] * z1wlo * z2wlo;

	result /= totweight;

	return result;

}
#endif // end brgastro::add_cache functions

const int brgastro::tfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = -brgastro::integrate_ltd( 0, brgastro::zfa( in_params ) ) / c;
	}
	catch(std::exception &e)
	{
		std::cerr << "ERROR: Could not calculate cache for " << _name_base() << "\n"
				<< "Exception: " << e.what() << "\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

// brgastro::tNFW_sig_cache class methods
#if (1)

const int brgastro::tNFW_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::tNFW_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		if ( make_array3d( signal, halo_z_res, halo_m_res, r_res ) )
			return 1;

		// Read in data
		i = j = k = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k] = temp_data;
			k++;
			if ( k >= r_res )
			{
				k = 0;
				j++;
				if ( j >= halo_m_res )
				{
					j = 0;
					i++;
					if ( i >= halo_z_res )
						break;
				}
			}
		}

		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_sig_cache::unload()
{
	loaded = false;
	return del_array3d( signal, halo_z_res, halo_m_res );
}

const int brgastro::tNFW_sig_cache::calc( const bool silent )
{
	int i, j, k;
	double mass, conc, tau;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	if ( make_array3d( signal, halo_z_res, halo_m_res, r_res ) )
		return 1;

	i = j = k = 0;

	for ( double z = halo_z_min; z < halo_z_max; z += halo_z_step )
	{
		for ( double m = halo_m_min; m < halo_m_max; m += halo_m_step )
		{
			mass = std::pow( 10, m ) * unitconv::Msuntokg;
			conc = brgastro::cfm( mass, z );
			tau = default_tau_factor * conc;
			for ( double r = r_min; r < r_max; r += r_step )
			{

				signal[i][j][k] =
						brgastro::tNFW_profile( mass, z, conc, tau ).deltasigma(
								r );
				k++;
			}
			k = 0;
			j++;
		}
		j = k = 0;
		i++;
	}

	return 0;
}

const int brgastro::tNFW_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				if ( !( out_file << signal[i][j][k] << "\t" ) )
					return errorNOS();
			}
			out_file << "\n";
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_sig_cache::set_file_name( const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		return unload();
	}
	return 0;
}

const int brgastro::tNFW_sig_cache::set_range( const double new_halo_z_min,
		const double new_halo_z_max, const double new_halo_z_step,
		const double new_halo_m_min, const double new_halo_m_max,
		const double new_halo_m_step, const double new_r_min,
		const double new_r_max, const double new_r_step, const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_sig_cache::set_precision( const int new_precision,
		const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr << "ERROR: Precision for tNFW_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const bool silent )
{
	double zlo, zhi, mlo, mhi, rlo, rhi;
	int z_i, m_i, r_i; // Lower nearby array points
	BRG_DISTANCE rwlo, rwhi;
	double mwlo, mwhi, zwlo, zwhi;
	double lm = log10( m * unitconv::kgtoMsun );
	double result = 0;

	if ( !loaded )
	{
		#pragma omp critical(tNFW_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}

		if ( result == -1 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: Could neither load nor calculate tNFW_sig_cache!\n";
			return result;
		}
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = lm - mlo;
	mwhi = mhi - lm;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	double totweight = halo_z_step * halo_m_step * r_step;

#ifdef _BRG_USE_UNITS_
	result = BRG_UNITS(0,-2,0,1,0,0,0); // To get the right units
#else
	result = 0;
#endif

	result += signal[z_i][m_i][r_i] * zwhi * mwhi * rwhi;
	result += signal[z_i][m_i][r_i + 1] * zwhi * mwhi * rwlo;
	result += signal[z_i][m_i + 1][r_i] * zwhi * mwlo * rwhi;
	result += signal[z_i][m_i + 1][r_i + 1] * zwhi * mwlo * rwlo;
	result += signal[z_i + 1][m_i][r_i] * zwlo * mwhi * rwhi;
	result += signal[z_i + 1][m_i][r_i + 1] * zwlo * mwhi * rwlo;
	result += signal[z_i + 1][m_i + 1][r_i] * zwlo * mwlo * rwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1] * zwlo * mwlo * rwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::tNFW_sig_cache functions

// brgastro::NFW_offset_sig_cache class methods
#if (1)

const int brgastro::tNFW_offset_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k, l;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::NFW_offset_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step >> offset_r_min >> offset_r_max >> offset_r_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		offset_r_res = (int)( ( offset_r_max - offset_r_min ) / offset_r_step )
				+ 1;
		if ( make_array4d( signal, halo_z_res, halo_m_res, r_res,
				offset_r_res ) )
			return 1;

		// Read in data
		i = j = k = l = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k][l] = temp_data;
			l++;
			if ( l >= offset_r_res )
			{
				l = 0;
				k++;
				if ( k >= r_res )
				{
					k = 0;
					j++;
					if ( j >= halo_m_res )
					{
						j = 0;
						i++;
						if ( i >= halo_z_res )
							break;
					}
				}
			}
		}
		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::unload()
{
	loaded = false;
	return del_array4d( signal, halo_z_res, halo_m_res, r_res );
}

const int brgastro::tNFW_offset_sig_cache::calc( const bool silent )
{
	double z, lm, r, l_offset_r;
	double mass, conc, tau, offset_r;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 )
			|| ( offset_r_max <= offset_r_min ) || ( offset_r_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_offset_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	offset_r_res = (int)( ( offset_r_max - offset_r_min ) / offset_r_step )
			+ 1;
	if ( make_array4d( signal, halo_z_res, halo_m_res, r_res, offset_r_res ) )
		return 1;

	z = halo_z_min;
	for ( int i = 0; i < halo_z_res; i++ )
	{
		lm = halo_m_min;
		for ( int j = 0; j < halo_m_res; j++ )
		{
			mass = std::pow( 10, lm ) * unitconv::Msuntokg;
			conc = brgastro::cfm( mass, z );
			tau = default_tau_factor * conc;
			r = r_min;
			for ( int k = 0; k < r_res; k++ )
			{
				l_offset_r = offset_r_min;
				for ( int l = 0; l < offset_r_res; l++ )
				{
					offset_r = std::pow( 10, l_offset_r ) * unitconv::kpctom;
					signal[i][j][k][l] = brgastro::tNFW_profile( mass, z, conc,
							tau ).offset_WLsig( r, offset_r );
					l_offset_r += offset_r_step;
				}
				r += r_step;
			}
			lm += halo_m_step;
		}
		z += halo_z_step;
	}

	return 0;
}

const int brgastro::tNFW_offset_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n" << offset_r_min
			<< "\t" << offset_r_max << "\t" << offset_r_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				for ( int l = 0; l < offset_r_res; l++ )
				{
					if ( !( out_file << signal[i][j][k][l] << "\n" ) )
						return errorNOS();
				}
			}
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_offset_sig_cache::set_file_name(
		const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		if ( unload() )
			return errorNOS();
	}
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::set_range(
		const double new_halo_z_min, const double new_halo_z_max,
		const double new_halo_z_step, const double new_halo_m_min,
		const double new_halo_m_max, const double new_halo_m_step,
		const double new_r_min, const double new_r_max,
		const double new_r_step, const double new_offset_r_min,
		const double new_offset_r_max, const double new_offset_r_step,
		const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step )
			|| ( halo_z_min != new_offset_r_min )
			|| ( halo_z_max != new_offset_r_max )
			|| ( halo_z_step != new_offset_r_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;
		offset_r_min = new_offset_r_min;
		offset_r_max = new_offset_r_max;
		offset_r_step = new_offset_r_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_offset_sig_cache::set_precision(
		const int new_precision, const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Precision for tNFW_offset_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_offset_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const BRG_DISTANCE offset_r,
		const bool silent )
{
	double zlo, zhi, mlo, mhi;
	BRG_DISTANCE rlo, rhi, orlo, orhi;
	int z_i, m_i, r_i, or_i; // Lower nearby array points
	double zwlo, zwhi, mwlo, mwhi, rwlo, rwhi, orwlo, orwhi;
	double lm = log10( m * unitconv::kgtoMsun );
	double lor = log10( offset_r * unitconv::mtokpc );
	double result = 0;

	if ( !loaded )
	{
		#pragma omp critical(tNFW_offset_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
			return result;
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = lm - mlo;
	mwhi = mhi - lm;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	or_i = (int)( ( lor - offset_r_min ) / offset_r_step );
	or_i = max( or_i, 0 );
	or_i = min( or_i, offset_r_res - 2 );

	orlo = offset_r_min + offset_r_step * r_i;
	orhi = offset_r_min + offset_r_step * ( r_i + 1 );
	orwlo = lor - orlo;
	orwhi = orhi - lor;

	double totweight = halo_z_step * halo_m_step * r_step * offset_r_step;

#ifdef _BRG_USE_UNITS_
	result = BRG_UNITS(0,-2,0,1,0,0,0);
#else
	result = 0;
#endif

	result += signal[z_i][m_i][r_i][or_i] * zwhi * mwhi * rwhi * orwhi;
	result += signal[z_i][m_i][r_i + 1][or_i] * zwhi * mwhi * rwlo * orwhi;
	result += signal[z_i][m_i + 1][r_i][or_i] * zwhi * mwlo * rwhi * orwhi;
	result += signal[z_i][m_i + 1][r_i + 1][or_i] * zwhi * mwlo * rwlo * orwhi;
	result += signal[z_i + 1][m_i][r_i][or_i] * zwlo * mwhi * rwhi * orwhi;
	result += signal[z_i + 1][m_i][r_i + 1][or_i] * zwlo * mwhi * rwlo * orwhi;
	result += signal[z_i + 1][m_i + 1][r_i][or_i] * zwlo * mwlo * rwhi * orwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1][or_i] * zwlo * mwlo * rwlo
			* orwhi;
	result += signal[z_i][m_i][r_i][or_i + 1] * zwhi * mwhi * rwhi * orwlo;
	result += signal[z_i][m_i][r_i + 1][or_i + 1] * zwhi * mwhi * rwlo * orwlo;
	result += signal[z_i][m_i + 1][r_i][or_i + 1] * zwhi * mwlo * rwhi * orwlo;
	result += signal[z_i][m_i + 1][r_i + 1][or_i + 1] * zwhi * mwlo * rwlo
			* orwlo;
	result += signal[z_i + 1][m_i][r_i][or_i + 1] * zwlo * mwhi * rwhi * orwlo;
	result += signal[z_i + 1][m_i][r_i + 1][or_i + 1] * zwlo * mwhi * rwlo
			* orwlo;
	result += signal[z_i + 1][m_i + 1][r_i][or_i + 1] * zwlo * mwlo * rwhi
			* orwlo;
	result += signal[z_i + 1][m_i + 1][r_i + 1][or_i + 1] * zwlo * mwlo * rwlo
			* orwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::NFW_offset_sig_cache functions

// brgastro::tNFW_group_sig_cache class methods
#if (1)

const int brgastro::tNFW_group_sig_cache::load( const bool silent )
{
	std::ifstream in_file;
	std::string file_data;
	bool need_to_calc = false;
	int loop_counter = 0;
	double temp_data;
	int i, j, k, l;

	if ( loaded )
		return 0;

	do
	{
		if ( loop_counter >= 2 )
		{
			if ( !silent )
				std::cerr
						<< "ERROR: infinite loop detected in brgastro::NFW_offset_sig_cache.\n";
			return INFINITE_LOOP_ERROR;
		}
		else
		{
			loop_counter++;
		}
		need_to_calc = false;

		if ( open_file( in_file, file_name ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Check that it has the right header
		getline( in_file, file_data );
		if ( file_data.compare( header_string ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Trim out any other commented lines
		if ( trim_comments_all_at_top( in_file ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Load range parameters;
		if ( !( in_file >> halo_z_min >> halo_z_max >> halo_z_step
				>> halo_m_min >> halo_m_max >> halo_m_step >> r_min >> r_max
				>> r_step >> group_c_min >> group_c_max >> group_c_step ) )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

		// Set up data
		halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
		halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_m_step ) + 1;
		r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
		group_c_res = (int)( ( group_c_max - group_c_min ) / group_c_step )
				+ 1;
		if ( make_array4d( signal, halo_z_res, halo_m_res, r_res,
				group_c_res ) )
			return 1;

		// Read in data
		i = j = k = l = 0;
		while ( !in_file.eof() )
		{
			in_file >> temp_data;
			signal[i][j][k][l] = temp_data;
			l++;
			if ( l >= group_c_res )
			{
				l = 0;
				k++;
				if ( k >= r_res )
				{
					k = 0;
					j++;
					if ( j >= halo_m_res )
					{
						j = 0;
						i++;
						if ( i >= halo_z_res )
							break;
					}
				}
			}
		}
		// Check that it was all read properly
		if ( i < halo_z_res )
		{
			need_to_calc = true;
			if ( calc( silent ) )
				return 1;
			unload();
			continue;
		}

	} while ( need_to_calc );

	// Finish up
	in_file.close();
	in_file.clear();
	loaded = true;
	return 0;
}

const int brgastro::tNFW_group_sig_cache::unload()
{
	loaded = false;
	return del_array4d( signal, halo_z_res, halo_m_res, r_res );
}

const int brgastro::tNFW_group_sig_cache::calc( const bool silent )
{
	double z, lm, r, l_group_c;
	double mass, group_c;

	// Test that range is sane
	if ( ( halo_z_max <= halo_z_min ) || ( halo_z_step <= 0 )
			|| ( halo_m_max <= halo_m_min ) || ( halo_m_step <= 0 )
			|| ( r_max <= r_min ) || ( r_step <= 0 )
			|| ( group_c_max <= group_c_min ) || ( group_c_step <= 0 ) )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Bad range passed to tNFW_group_sig_cache::calc(silent)\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	// Set up data
	halo_z_res = (int)( ( halo_z_max - halo_z_min ) / halo_z_step ) + 1;
	halo_m_res = (int)( ( halo_m_max - halo_m_min ) / halo_z_step ) + 1;
	r_res = (int)( ( r_max - r_min ) / r_step ) + 1;
	group_c_res = (int)( ( group_c_max - group_c_min ) / group_c_step ) + 1;
	if ( make_array4d( signal, halo_z_res, halo_m_res, r_res, group_c_res ) )
		return 1;

	z = halo_z_min;
	for ( int i = 0; i < halo_z_res; i++ )
	{
		lm = halo_m_min;
		for ( int j = 0; j < halo_m_res; j++ )
		{
			mass = std::pow( 10, lm ) * unitconv::Msuntokg;
			r = r_min;
			for ( int k = 0; k < r_res; k++ )
			{
				l_group_c = group_c_min;
				for ( int l = 0; l < group_c_res; l++ )
				{
					group_c = std::pow( 10, l_group_c ) * unitconv::kpctom;
					signal[i][j][k][l] =
							brgastro::tNFW_profile( mass, z ).group_WLsig( r,
									group_c );
					l_group_c += group_c_step;
				}
				r += r_step;
			}
			lm += halo_m_step;
		}
		z += halo_z_step;
	}

	return 0;
}

const int brgastro::tNFW_group_sig_cache::output( const bool silent )
{
	std::ofstream out_file;
	std::string file_data;

	if ( !loaded )
	{
		if ( calc( silent ) )
			return 1;
	}

	if ( open_file( out_file, file_name ) )
		return 1;

	// Output header
	out_file << header_string << "\n#\n";

	// Output range
	out_file << halo_z_min << "\t" << halo_z_max << "\t" << halo_z_step << "\n"
			<< halo_m_min << "\t" << halo_m_max << "\t" << halo_m_step << "\n"
			<< r_min << "\t" << r_max << "\t" << r_step << "\n" << group_c_min
			<< "\t" << group_c_max << "\t" << group_c_step << "\n";

	// Output data
	for ( int i = 0; i < halo_z_res; i++ )
	{
		for ( int j = 0; j < halo_m_res; j++ )
		{
			for ( int k = 0; k < r_res; k++ )
			{
				for ( int l = 0; l < group_c_res; l++ )
				{
					if ( !( out_file << signal[i][j][k][l] << "\n" ) )
						return FILE_ACCESS_ERROR;
				}
			}
		}
	}

	out_file.close();
	out_file.clear();

	return 0;

}

const int brgastro::tNFW_group_sig_cache::set_file_name(
		const std::string new_name )
{
	file_name = new_name;
	if ( loaded )
	{
		if ( unload() )
			return errorNOS();
	}
	return 0;
}

const int brgastro::tNFW_group_sig_cache::set_range(
		const double new_halo_z_min, const double new_halo_z_max,
		const double new_halo_z_step, const double new_halo_m_min,
		const double new_halo_m_max, const double new_halo_m_step,
		const double new_r_min, const double new_r_max,
		const double new_r_step, const double new_group_c_min,
		const double new_group_c_max, const double new_group_c_step,
		const bool silent )
{
	// First we try to load
	if ( !loaded )
		load( silent );

	// Go through variables, check if any are actually changed. If so, recalculate cache
	if ( ( halo_z_min != new_halo_z_min ) || ( halo_z_max != new_halo_z_max )
			|| ( halo_z_step != new_halo_z_step )
			|| ( halo_m_min != new_halo_m_min )
			|| ( halo_m_max != new_halo_m_max )
			|| ( halo_m_step != new_halo_m_step )
			|| ( halo_z_min != new_r_min ) || ( halo_z_max != new_r_max )
			|| ( halo_z_step != new_r_step )
			|| ( halo_z_min != new_group_c_min )
			|| ( halo_z_max != new_group_c_max )
			|| ( halo_z_step != new_group_c_step ) )
	{
		halo_z_min = new_halo_z_min;
		halo_z_max = new_halo_z_max;
		halo_z_step = new_halo_z_step;
		halo_m_min = new_halo_m_min;
		halo_m_max = new_halo_m_max;
		halo_m_step = new_halo_m_step;
		r_min = new_r_min;
		r_max = new_r_max;
		r_step = new_r_step;
		group_c_min = new_group_c_min;
		group_c_max = new_group_c_max;
		group_c_step = new_group_c_step;

		if ( unload() )
			return errorNOS( silent );
		if ( calc( silent ) )
			return 1;
	}
	return 0;
}

const int brgastro::tNFW_group_sig_cache::set_precision(
		const int new_precision, const bool silent )
{
	if ( new_precision > 0 )
	{
		sig_digits = min( new_precision, DBL_MAX_PRECISION );
		return 0;
	}
	else
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Precision for tNFW_group_sig_cache must be > 0.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
}

const BRG_UNITS brgastro::tNFW_group_sig_cache::get( const double z,
		const BRG_MASS m, const BRG_DISTANCE r, const double group_c,
		const bool silent )
{
	double zlo, zhi, mlo, mhi, rlo, rhi, gclo, gchi;
	int z_i, m_i, r_i, gc_i; // Lower nearby array points
	double zwlo, zwhi, mwlo, mwhi, rwlo, rwhi, gcwlo, gcwhi;
	double lm = log10( m * unitconv::kgtoMsun );
#ifdef _BRG_USE_UNITS_
	BRG_UNITS result(0,-2,0,1,0,0,0);
#else
	double result = 0;
#endif

	if ( !loaded )
	{
		#pragma omp critical(tNFW_group_sig_load)
		{
			if ( load( silent ) )
			{
				result = -1;
			}
		}
		if ( result == -1 )
			return result;
	}

	z_i = (int)( ( z - halo_z_min ) / halo_z_step );
	z_i = max( z_i, 0 );
	z_i = min( z_i, halo_z_res - 2 );

	zlo = halo_z_min + halo_z_step * z_i;
	zhi = halo_z_min + halo_z_step * ( z_i + 1 );
	zwlo = z - zlo;
	zwhi = zhi - z;

	m_i = (int)( ( lm - halo_m_min ) / halo_m_step );
	m_i = max( m_i, 0 );
	m_i = min( m_i, halo_m_res - 2 );

	mlo = halo_m_min + halo_m_step * m_i;
	mhi = halo_m_min + halo_m_step * ( m_i + 1 );
	mwlo = m - mlo;
	mwhi = mhi - m;

	r_i = (int)( ( r - r_min ) / r_step );
	r_i = max( r_i, 0 );
	r_i = min( r_i, r_res - 2 );

	rlo = r_min + r_step * r_i;
	rhi = r_min + r_step * ( r_i + 1 );
	rwlo = r - rlo;
	rwhi = rhi - r;

	gc_i = (int)( ( group_c - group_c_min ) / group_c_step );
	gc_i = max( gc_i, 0 );
	gc_i = min( gc_i, group_c_res - 2 );

	gclo = group_c_min + group_c_step * r_i;
	gchi = group_c_min + group_c_step * ( r_i + 1 );
	gcwlo = group_c - gclo;
	gcwhi = gchi - group_c;

	double totweight = halo_z_step * halo_m_step * r_step * group_c_step;

	result = 0;

	result += signal[z_i][m_i][r_i][gc_i] * zwhi * mwhi * rwhi * gcwhi;
	result += signal[z_i][m_i][r_i + 1][gc_i] * zwhi * mwhi * rwlo * gcwhi;
	result += signal[z_i][m_i + 1][r_i][gc_i] * zwhi * mwlo * rwhi * gcwhi;
	result += signal[z_i][m_i + 1][r_i + 1][gc_i] * zwhi * mwlo * rwlo * gcwhi;
	result += signal[z_i + 1][m_i][r_i][gc_i] * zwlo * mwhi * rwhi * gcwhi;
	result += signal[z_i + 1][m_i][r_i + 1][gc_i] * zwlo * mwhi * rwlo * gcwhi;
	result += signal[z_i + 1][m_i + 1][r_i][gc_i] * zwlo * mwlo * rwhi * gcwhi;
	result += signal[z_i + 1][m_i + 1][r_i + 1][gc_i] * zwlo * mwlo * rwlo
			* gcwhi;
	result += signal[z_i][m_i][r_i][gc_i + 1] * zwhi * mwhi * rwhi * gcwlo;
	result += signal[z_i][m_i][r_i + 1][gc_i + 1] * zwhi * mwhi * rwlo * gcwlo;
	result += signal[z_i][m_i + 1][r_i][gc_i + 1] * zwhi * mwlo * rwhi * gcwlo;
	result += signal[z_i][m_i + 1][r_i + 1][gc_i + 1] * zwhi * mwlo * rwlo
			* gcwlo;
	result += signal[z_i + 1][m_i][r_i][gc_i + 1] * zwlo * mwhi * rwhi * gcwlo;
	result += signal[z_i + 1][m_i][r_i + 1][gc_i + 1] * zwlo * mwhi * rwlo
			* gcwlo;
	result += signal[z_i + 1][m_i + 1][r_i][gc_i + 1] * zwlo * mwlo * rwhi
			* gcwlo;
	result += signal[z_i + 1][m_i + 1][r_i + 1][gc_i + 1] * zwlo * mwlo * rwlo
			* gcwlo;

	result /= totweight;

	return result;

}

#endif // end brgastro::tNFW_group_sig_cache functions

// brgastro::accel_function class methods
#if (1)

const int brgastro::accel_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::accel_function::operator()( const BRG_UNITS & in_param,
BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to accel_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = _host_ptr_->accel( in_param, silent );
	return 0;
}

brgastro::accel_function::accel_function()
{
	_host_ptr_ = 0;
}
brgastro::accel_function::accel_function( const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::accel_function function implementations

// brgastro::spherical_density_function class methods
#if (1)

const int brgastro::spherical_density_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::spherical_density_function::operator()(
		const BRG_UNITS & in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to spherical_density_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = 4 * pi * std::pow( in_param, 2 )
			* _host_ptr_->dens( in_param );
	return 0;
}

brgastro::spherical_density_function::spherical_density_function()
{
	_host_ptr_ = 0;
}
brgastro::spherical_density_function::spherical_density_function(
		const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::spherical_density_function class methods

// brgastro::projected_density_function class methods
#if (1)

const int brgastro::projected_density_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::projected_density_function::set_offset_R(
		const BRG_DISTANCE &new_offset_R )
{
	_offset_R_ = new_offset_R;
	return 0;
}

const int brgastro::projected_density_function::operator()(
		const BRG_UNITS & in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to projected_density_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	BRG_DISTANCE r = quad_add( in_param, _offset_R_ );
	out_param = _host_ptr_->dens( r );
	return 0;
}

brgastro::projected_density_function::projected_density_function()
{
	_host_ptr_ = 0;
	_offset_R_ = 0;
}

brgastro::projected_density_function::projected_density_function(
		const density_profile *init_host, const BRG_DISTANCE &init_offset_R )
{
	set_host_ptr( init_host );
	set_offset_R( init_offset_R );
}

#endif // end brgastro::projected_density_function class methods

// brgastro::cylindrical_density_function class methods
#if (1)

const int brgastro::cylindrical_density_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::cylindrical_density_function::operator()(
		const BRG_UNITS & in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to cylindrical_density_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	return 0;
}

brgastro::cylindrical_density_function::cylindrical_density_function()
{
	_host_ptr_ = 0;
}
brgastro::cylindrical_density_function::cylindrical_density_function(
		const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end brgastro::cylindrical_density_function class methods

// brgastro::solve_rhm_function class methods
#if (1)
const int brgastro::solve_rhm_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

const int brgastro::solve_rhm_function::set_target_mass(
		const BRG_MASS &new_target_mass )
{
	_target_mass_ = new_target_mass;
	return 0;
}

const int brgastro::solve_rhm_function::operator()( const BRG_UNITS & in_param,
BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to solve_rhm_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	double result = fabs(
			_target_mass_ - _host_ptr_->enc_mass( fabs( in_param ), silent ) );
	out_param = result;
	return 0;
}

brgastro::solve_rhm_function::solve_rhm_function()
{
	_host_ptr_ = 0;
	_target_mass_ = 0;
}
;
brgastro::solve_rhm_function::solve_rhm_function(
		const density_profile *init_host, const BRG_MASS &new_target_mass )
{
	set_host_ptr( init_host );
	set_target_mass( new_target_mass );
}
#endif // end brgastro::solve_rhm_function function implementations

// brgastro::offset_ring_dens_function class methods
#if (1)

const int brgastro::offset_ring_dens_function::set_R0(
		const BRG_DISTANCE &new_R0 )
{
	_R0_ = new_R0;
	return 0;
}

const int brgastro::offset_ring_dens_function::set_R(
		const BRG_DISTANCE &new_R )
{
	_R_ = new_R;
	return 0;
}

const int brgastro::offset_ring_dens_function::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_ring_dens_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	BRG_DISTANCE d = lc_add( _R0_, _R_, in_param );

	out_param = _host_ptr_->proj_dens( d, silent );

	return 0;
}

const int brgastro::offset_ring_dens_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_ring_dens_function::offset_ring_dens_function()
{
	_host_ptr_ = 0;
	_R_ = 0;
	_R0_ = 0;
}

brgastro::offset_ring_dens_function::offset_ring_dens_function(
		const density_profile *new_host, const BRG_DISTANCE &new_R_0,
		const BRG_DISTANCE &new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
	_R0_ = new_R_0;
}

#endif

// brgastro::offset_circ_dens_function class methods
#if (1)

const int brgastro::offset_circ_dens_function::set_R0(
		const BRG_DISTANCE &new_R_0 )
{
	_R0_ = new_R_0;
	return 0;
}

const int brgastro::offset_circ_dens_function::operator()(
		const std::vector< BRG_UNITS > &in_params,
		std::vector< BRG_UNITS > & out_params, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_circ_dens_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	if ( in_params.size() != 2 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: in_params.size() must == 2 in offset_circ_dens_function.\n";
		return INVALID_ARGUMENTS_ERROR;
	}
	out_params.clear();

	BRG_DISTANCE R = in_params[0];
	BRG_DISTANCE d = lc_add( _R0_, R, in_params[1] );

	BRG_UNITS result = _host_ptr_->proj_dens( d, silent );
	out_params.push_back( result );

	return 0;
}

const int brgastro::offset_circ_dens_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_circ_dens_function::offset_circ_dens_function()
{
	_host_ptr_ = 0;
	_R0_ = 0;
}

brgastro::offset_circ_dens_function::offset_circ_dens_function(
		const density_profile *new_host, const BRG_DISTANCE &new_R_0 )
{
	_host_ptr_ = new_host;
	_R0_ = new_R_0;
}
#endif // end brgastro::offset_circ_dens_function function implementations

// brgastro::offset_WLsig_function class methods
#if (1)

const int brgastro::offset_WLsig_function::set_R( const BRG_DISTANCE &new_R )
{
	_R_ = new_R;
	return 0;
}

const int brgastro::offset_WLsig_function::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_WLsig_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	out_param = _host_ptr_->offset_WLsig( _R_, in_param, silent );

	return 0;
}

const int brgastro::offset_WLsig_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_WLsig_function::offset_WLsig_function()
{
	_host_ptr_ = 0;
	_R_ = 0;
}

brgastro::offset_WLsig_function::offset_WLsig_function(
		const density_profile *new_host, const BRG_DISTANCE &new_R )
{
	_host_ptr_ = new_host;
	_R_ = new_R;
}

#endif

// brgastro::offset_WLsig_weight_function class methods
#if (1)

const int brgastro::offset_WLsig_weight_function::set_c( const double new_c )
{
	_c_ = new_c;
	return 0;
}

const int brgastro::offset_WLsig_weight_function::operator()(
		const BRG_UNITS &in_param,
		BRG_UNITS & out_param, const bool silent ) const
{
	BRG_UNITS result;
	if ( _host_ptr_ == 0 )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to offset_WLsig_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}

	if ( _c_ == 0 )
	{
		result = 2 * pi * in_param * _host_ptr_->proj_dens( in_param );
	}
	else
	{
		BRG_UNIQUE_PTR<brgastro::density_profile> group_profile(_host_ptr_->density_profile_clone());
		group_profile->set_c(_c_);
		result = 2 * pi * in_param * group_profile->proj_dens(in_param);
		del_obj(group_profile);
	}

	out_param = result;

	return 0;
}

const int brgastro::offset_WLsig_weight_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

brgastro::offset_WLsig_weight_function::offset_WLsig_weight_function()
{
	_host_ptr_ = 0;
	_c_ = 0;
}

brgastro::offset_WLsig_weight_function::offset_WLsig_weight_function(
		const density_profile *new_host, const double init_c )
{
	_host_ptr_ = new_host;
	_c_ = init_c;
}

#endif

#endif// end Class methods

/** Global Function Definitions **/
#if (1)
// grid functions
#if (1)
const int brgastro::get_ra_grid( const BRG_ANGLE &ra )
{
	grid_cache cache;
	return (int)floor(
			( ra - cache.ra_grid_min() ) / safe_d( cache.ra_grid_step() ) );
}
const int brgastro::get_dec_grid( const BRG_ANGLE &dec )
{
	grid_cache cache;
	return (int)floor( ( dec - cache.dec_grid_min() ) / cache.dec_grid_step() );
}
const int brgastro::get_z_grid( const double z )
{
	grid_cache cache;
	return (int)floor( ( z - cache.z_grid_min() ) / cache.z_grid_step() );
}

const BRG_ANGLE brgastro::get_ra_grid_lower( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ra_grid;
}
const BRG_ANGLE brgastro::get_dec_grid_lower( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * dec_grid;
}
const double brgastro::get_z_grid_lower( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * z_grid;
}

const BRG_ANGLE brgastro::get_ra_grid_upper( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ( ra_grid + 1 );
}
const BRG_ANGLE brgastro::get_dec_grid_upper( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * ( dec_grid + 1 );
}
const double brgastro::get_z_grid_upper( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * ( z_grid + 1 );
}

const BRG_ANGLE brgastro::get_ra_grid_mid( const int ra_grid )
{
	grid_cache cache;
	return cache.ra_grid_min() + cache.ra_grid_step() * ( ra_grid + 0.5 );
}
const BRG_ANGLE brgastro::get_dec_grid_mid( const int dec_grid )
{
	grid_cache cache;
	return cache.dec_grid_min() + cache.dec_grid_step() * ( dec_grid + 0.5 );
}
const double brgastro::get_z_grid_mid( const int z_grid )
{
	grid_cache cache;
	return cache.z_grid_min() + cache.z_grid_step() * ( z_grid + 0.5 );
}
#endif // end grid functions

// dfa and afd functions
#if (1)

const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &da, const double z )
{
	return da * dfa_cache().get( z );
}
const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &a1, const BRG_ANGLE &a2,
		const double z )
{
	return brgastro::dfa( a2 - a1, z );
}
const BRG_DISTANCE brgastro::dfa( const BRG_ANGLE &a1x, const BRG_ANGLE &a1y,
		const BRG_ANGLE &a2x, const BRG_ANGLE &a2y, const double z )
{
	return brgastro::dfa( skydist2d( a1x, a1y, a2x, a2y ), z );
}
const BRG_DISTANCE brgastro::dfa( const brgastro::sky_obj *obj1,
		const brgastro::sky_obj *obj2, const double z )
{
	double z_to_use;
	if ( z == -1 )
		z_to_use = obj1->z();
	else
		z_to_use = z;
	return brgastro::dfa(
			skydist2d( obj1->ra(), obj1->dec(), obj2->ra(), obj2->dec() ),
			z_to_use );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &dd, const double z )
{
	return dd / safe_d( dfa_cache().get( z ) );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &d1, const BRG_DISTANCE &d2,
		const double z )
{
	return brgastro::afd( fabs( d2 - d1 ), z );
}
const BRG_ANGLE brgastro::afd( const BRG_DISTANCE &d1x,
		const BRG_DISTANCE &d1y, const BRG_DISTANCE &d2x,
		const BRG_DISTANCE &d2y, const double z )
{
	return brgastro::afd( dist2d( d1x, d1y, d2x, d2y ), z );
}

const double brgastro::zfa( const double a )
{
	return 1. / safe_d( a ) - 1.;
}
const double brgastro::afz( const double z )
{
	return 1. / safe_d( 1 + z );
}

const BRG_TIME brgastro::tfz( const double z )
{
	return brgastro::tfa( afz( z ) );
}
const BRG_TIME brgastro::tfa( const double a )
{
	return tfa_cache().get( a );
}
const double brgastro::zft( const BRG_TIME &t )
{
	return brgastro::zfa( brgastro::aft( t ) );
}
const double brgastro::aft( const BRG_TIME &t )
{
	brgastro::tfa_cache cache;
	return cache.inverse_get( t );
}

#endif

// Functions to integrate out distances
#if (1)
const double brgastro::integrate_add( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 0 );
}
const double brgastro::integrate_cmd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 1 );
}
const double brgastro::integrate_Ld( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 2 );
}
const double brgastro::integrate_ltd( const double z1, const double z2 )
{
	return brgastro::integrate_distance( z1, z2, 3, 10000000 );
}
const double brgastro::integrate_add( const double z )
{
	return brgastro::integrate_distance( 0, z, 0 );
}
const double brgastro::integrate_cmd( const double z )
{
	return brgastro::integrate_distance( 0, z, 1 );
}
const double brgastro::integrate_Ld( const double z )
{
	return brgastro::integrate_distance( 0, z, 2 );
}
const double brgastro::integrate_ltd( const double z )
{
	return brgastro::integrate_distance( 0, z, 3 );
}
const double brgastro::integrate_distance( const double z1_init,
		const double z2_init, const int mode, const int n )
{
	// Function that will help calculate cosmic distances, thanks to Richard Powell - http://www.atlasoftheuniverse.com/
	// NOTE: Will return a negative value if z2 < z1. This is by design.

	double OK0;
	double OM0 = Omega_m, OR0 = Omega_r, OL0 = Omega_l;
	double OM, OR, OL, OK;
	double z1 = z1_init, z2 = z2_init;
	double HD; //Hubble distance in billions of lightyears
	double z, a, a1, a2, adot, h1;
	double DC = DBL_MAX, DCC = 0, DT = DBL_MAX, DTT = 0, DM;
	//double age, size;
	int i;
	short int sign = 1;

	OK0 = 1 - OM0 - OL0 - OR0;

	if ( z2 < z1 )
	{
		std::swap(z1,z2);
		sign = -1;
	}

	if ( z1 == 0 )
	{
		z = z2;
		h1 = H_0;
		OM = OM0;
		OR = OR0;
		OL = OL0;
		OK = OK0;
	}
	else
	{
		a1 = 1 / ( 1 + z1 );
		a2 = 1 / ( 1 + z2 );
		z = ( a1 / a2 ) - 1;
		h1 = H_0
				* sqrt(
						OR0 * std::pow( a1, -4 ) + OM0 * std::pow( a1, -3 )
								+ OK0 * std::pow( a1, -2 ) + OL0 );
		OM = OM0 * std::pow( H_0 / h1, 2 ) * std::pow( a1, -3 );
		OR = OR0 * std::pow( H_0 / h1, 2 ) * std::pow( a1, -4 );
		OL = OL0 * std::pow( H_0 / h1, 2 );
		OK = 1 - OM - OR - OL;
	}

	HD = c / h1;

	for ( i = n; i >= 1; i-- )        // This loop is the numerical integration
	{
		a = ( i - 0.5 ) / n;              // Steadily decrease the scale factor
		// Comoving formula (See section 4 of Hogg, but I've added a radiation term too):
		adot = a
				* sqrt(
						OM * std::pow( 1 / a, 3 ) + OK * std::pow( 1 / a, 2 )
								+ OL + OR * std::pow( 1 / a, 4 ) ); // Note that "a" is equivalent to 1/(1+z)
		 // Collect DC and DT until the correct scale is reached
		DCC = DCC + 1 / ( a * adot ) / n; // Running total of the comoving distance
		DTT = DTT + 1 / adot / n; // Running total of the light travel time (see section 10 of Hogg)
		 // Collect DC and DT until the correct scale is reached
		DC = DCC;                 // Comoving distance DC
		DT = DTT;                 // Light travel time DT
		if ( a <= 1 / ( 1 + z ) ) // Collect DC and DT until the correct scale is reached
		{
			break;
		}
	}

	// Transverse comoving distance DM from section 5 of Hogg:
	if ( OK > 0.0001 )
		DM = ( 1 / sqrt( OK ) ) * sinh( sqrt( OK ) * DC );
	else if ( OK < -0.0001 )
		DM = ( 1 / sqrt( fabs( OK ) ) ) * sin( sqrt( fabs( OK ) ) * DC );
	else
		DM = DC;

	//age = HD*DTT;                 // Age of the universe (billions of years)
	//size = HD*DCC;                // Comoving radius of the observable universe

	switch ( mode )
	{
	case 0:
		return sign * HD * DM / ( 1 + z ); // Angular diameter distance (section 6 of Hogg)
		break;
	case 1:
		return sign * HD * DC;             // Comoving distance
		break;
	case 2:
		return sign * HD * DM * ( 1 + z );  // Luminosity distance (section 7 of Hogg)
		break;
	case 3:
		return sign * HD * DT;             // Light travel distance
		break;
	default:
		return sign * HD * DT;             // Light travel distance
	}
}
#endif

// Other functions
#if (1)

const double brgastro::taufm( const double m_ratio, double conc,
		double tau_init, const bool silent )
{

	//Gives tau for a given Mtot/Mvir.
	if ( conc <= 0 )
	{
		conc = default_c;
		if ( !silent )
			std::cerr
					<< "WARNING: Invalid c value passed to taufm function.\n";
	}
	if ( tau_init <= 0 )
		tau_init = default_tau_factor * conc;
	double m_target = m_ratio * mftau( tau_init, conc );
	double taustepsize = tau_init / 2;
	double tautest[3];
	double mtest, mbest;
	int i, ibest;

	tautest[0] = tau_init / 2;
	mbest = 1e99;

	while ( taustepsize > 0.00001 * conc )
	{
		taustepsize /= 2;
		tautest[1] = tautest[0] - taustepsize;
		tautest[2] = tautest[0] + taustepsize;
		ibest = 0;
		for ( i = 0; i < 3; i++ )
		{
			mtest = mftau( tautest[i], conc );
			if ( fabs( mtest - m_target ) <= fabs( mbest - m_target ) )
			{
				ibest = i;
				mbest = mtest;
			}
		}
		tautest[0] = tautest[ibest];
	}

	return tautest[0];

}

const BRG_TIME brgastro::period( const brgastro::density_profile *host,
		const BRG_DISTANCE &r, const BRG_VELOCITY &vr, const BRG_VELOCITY &vt )
{
	BRG_UNITS mu = host->enc_mass( r ) * Gc;
	BRG_VELOCITY v = quad_add( vr, vt );
	BRG_DISTANCE a = -mu / 2 / safe_d( v * v / 2 - mu / safe_d( r ) );
	BRG_TIME result = (
			a > 0 ? 2 * pi * sqrt( std::pow( a, 3 ) / mu ) : BRG_TIME( 0 ) );
	return result;
}

const BRG_DISTANCE brgastro::ad_distance( double z1, double z2 )
{
	if ( z2 < z1 )
		std::swap( z1, z2 );
	return brgastro::add_cache().get( z1, z2 );
}

const BRG_UNITS brgastro::sigma_crit( const double z_lens,
		const double z_source )
{
	return pow( c, 2 ) / ( 4. * pi * Gc )
			* brgastro::ad_distance( 0, z_source )
			/ ( brgastro::ad_distance( 0, z_lens )
					* brgastro::ad_distance( z_lens, z_source ) );

}
#endif // end other functions

#endif // end Global function definitions
