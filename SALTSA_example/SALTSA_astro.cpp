#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>

#include "SALTSA_global.h"

#include "SALTSA_astro.h"
#include "SALTSA_misc_functions.hpp"
#include "SALTSA_calculus.hpp"
#include "SALTSA_solvers.hpp"

/** Static Class Initialisation **/
#if (1)

// Initialisation for SALTSA::tfa_cache
DEFINE_BRG_CACHE_STATIC_VARS( tfa_cache, 0.001, 1.02, 0.001 );

#endif // end Static Class Initialisation

/** Class Method Definitions **/
#if (1)
// SALTSA::redshift_obj class methods
#if (1)

const double SALTSA::redshift_obj::H() const
{
	// If not cached, calculate and cache it
	if(!_H_cached_)
	{
		if(_z_==0)
		{
			_H_cache_ = H_0;
		}
		else
		{
			_H_cache_ = SALTSA::H(_z_);
		}
		_H_cached_ = true;
	}
	return _H_cache_;
}

#endif

// SALTSA::density_profile class methods
#if (1)
/**
 *
 * @param r
 * @param silent
 * @return
 */
const double SALTSA::density_profile::Daccel( const double r,
		const bool silent ) const
{
	double dr;
	double a1, a2;
	// It's simplest, and a ton faster to just manually differentiate here.
	dr = max( r * SMALL_FACTOR, SMALL_FACTOR );
	a1 = accel( r, silent );
	a2 = accel( r + dr, silent );
	return ( a2 - a1 ) / safe_d( dr );
}

/**
 *
 * @param silent
 * @return
 */
const double SALTSA::density_profile::rhmtot( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmtot_cached )
		return _rhmtot_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmtot();
	solve_rhm_function func( this, target_mass );
	solve_rhm_function *func_ptr = &func;

	double max_r( default_tau_factor * rvir() );
	double rhm_test( 0 );

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

/**
 *
 * @param silent
 * @return
 */
const double SALTSA::density_profile::rhmvir( const bool silent ) const
{
	// If cached, return the cached value
	if ( hmvir_cached )
		return _rhmvir_cache_;

	// Not cached, so we'll have to calculate it
	double target_mass = hmvir();
	solve_rhm_function func( this, target_mass );
	solve_rhm_function *func_ptr = &func;

	double max_r( default_tau_factor * rvir() );
	double rhm_test( 0 );

	// First check for zero mass/radius/density
	if ( ( mvir() <= 0 ) || ( rvir() <= 0 ) || ( dens( rvir() / 2 ) < 0 ) )
	{
		hmvir_cached = true;
		return _rhmvir_cache_ = 0;
	}

	if ( solve_grid( func_ptr, 1, 0., max_r, 10, 0., rhm_test ) ) // If we can't solve it
	{
		if ( !silent )
			std::cerr << "WARNING: Could not solve half-mass radius. Assuming it's zero.\n";
		return -1;
	}
	else
	{
		_rhmvir_cache_ = max(0,rhm_test);
		hmvir_cached = true;
	}
	return _rhmvir_cache_;
}

/**
 *
 * @param r
 * @param silent
 * @return
 */
const double SALTSA::density_profile::enc_mass( const double r,
		const bool silent ) const
{
	if ( r == 0 )
		return 0;
	double r_to_use = std::fabs( r );
	SALTSA::spherical_density_function func( this );
	unsigned int num_in_params = 1, num_out_params = 0;
	double min_in_params( 0 ), max_in_params( r_to_use ), out_params( 0 );
	if ( SALTSA::integrate_Rhomberg( &func, num_in_params, min_in_params,
			max_in_params, num_out_params, out_params ) )
	{
		if ( !silent )
			std::cerr
					<< "WARNING: Could not integrate enclosed mass of density profile.\n";
	}
	return out_params;
}

#endif

// SALTSA::tNFW_profile class methods
#if (1)

const double SALTSA::tNFW_profile::_taufm( const double m_ratio,
		double precision, const bool silent ) const
{
	double m_target = m_ratio * _mftau();
	double taustepsize = _tau_ / 2;
	double tautest[3];
	double mtest, mbest;
	int i, ibest;

	tautest[0] = _tau_ / 2;
	mbest = 1e99;

	while ( taustepsize > precision * _c_ )
	{
		taustepsize /= 2;
		tautest[1] = tautest[0] - taustepsize;
		tautest[2] = tautest[0] + taustepsize;
		ibest = 0;
		for ( i = 0; i < 3; i++ )
		{
			mtest = _mftau( tautest[i] );
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

#if (1) // Constructors

/**
 *
 */
SALTSA::tNFW_profile::tNFW_profile()
{
	_mvir0_ = 0;
	_c_ = 0;
	_tau_ = 0;
}

/**
 *
 * @param init_mvir0
 * @param init_z
 * @param init_c
 * @param init_tau
 */
SALTSA::tNFW_profile::tNFW_profile( const double init_mvir0,
		const double init_z, const double init_c, const double init_tau ) :
		redshift_obj( init_z )
{
	_mvir0_ = init_mvir0;
	if ( init_c <= 0 )
	{
		_c_ = _cfm();
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
/**
 *
 */
SALTSA::tNFW_profile::~tNFW_profile()
{
}

#if (1) // Set functions

/**
 *
 * @param new_halo_mass
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::set_mvir( const double new_halo_mass,
		const bool silent )
{
	_mvir0_ = new_halo_mass;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
/**
 *
 * @param new_halo_tau
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::set_tau( const double new_halo_tau,
		const bool silent )
{
	_tau_ = new_halo_tau;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
/**
 *
 * @param new_halo_c
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::set_c( const double new_halo_c,
		const bool silent )
{
	_c_ = new_halo_c;
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
/**
 *
 * @param new_z
 * @return
 */
const int SALTSA::tNFW_profile::set_z( const double new_z )
{
	redshift_obj::set_z( new_z );
	hmvir_cached = false;
	hmtot_cached = false;
	return 0;
}
/**
 *
 * @param num_parameters
 * @param parameters
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< double > &parameters, const bool silent )
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
		if ( set_c( _cfm( parameters.at( 0 ), parameters.at( 1 ) ) ) )
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

/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::mvir() const
{
	return enc_mass(rvir()); // Not technically correct, but close enough for our purposes
}
/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::mvir0() const
{
	return _mvir0_;
}

/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::tau() const
{
	return _tau_;
}
/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::c() const
{
	return _c_;
}

/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::mtot() const
{
	return _mvir0_ * _mftau();
}

/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * _mvir0_, 1. / 3. );
}
/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::rvir() const
{
	return vvir() / H() / 10;
}
/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::rs() const
{
	return rvir() / _c_;
}
/**
 *
 * @param silent
 * @return
 */
const double SALTSA::tNFW_profile::rt( const bool silent ) const
{
	return rvir() / _tau_;
}

/**
 *
 * @return
 */
const double SALTSA::tNFW_profile::hmvir() const
{
	return enc_mass( rvir() ) / 2;
}


/**
 *
 * @param r
 * @return
 */
const double SALTSA::tNFW_profile::dens( const double r ) const
{
	double result, rho_c;

	double d_c, x, tau_use;
	if ( _tau_ <= 0 )
		tau_use = default_tau_factor * _c_;
	else
		tau_use = _tau_;
	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();

	result = ( d_c * rho_c ) / ( x * square( 1 + x ) )
			* square( tau_use )
			/ ( square( tau_use ) + square( x ) );

	return result;
}


/**
 *
 * @param r
 * @param silent
 * @return
 */
const double SALTSA::tNFW_profile::enc_mass( const double r,
		const bool silent ) const
{
	using SALTSA::square;
	using SALTSA::cube;
	using std::log;
	using std::atan;

	double rho_c;
	double m0, mx;
	double d_c, x, tau_use;
	if ( _tau_ < 0 )
		tau_use = default_tau_factor * _c_;
	else if (_tau_ == 0)
		return 0;
	else
		tau_use = _tau_;

	double tau_sq = square(tau_use);

	d_c = _delta_c();
	rho_c = 3 * square(H()) / ( 8 * pi * Gc );
	x = r / rs();

	// Result here integrated with Wolfram Alpha
	m0 = (2 * (1 + tau_sq) - (-1 + tau_sq) * 2 * log(tau_use));
	mx = ((2 * (1 + tau_sq)) / (1 + x)
							+ 4 * tau_use * atan(x / tau_use)
							+ 2 * (-1 + tau_sq) * log(1 + x)
							- (-1 + tau_sq) * log(tau_sq + square(x))) - m0;
	return cube(rs()) * d_c * rho_c * 2 * pi * tau_sq * mx
			/ square(1 + tau_sq);
}


/**
 *
 * @param parameters
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::get_parameters( std::vector< double > & parameters,
		const bool silent ) const
{
	parameters.resize( num_parameters() );

	try
	{
		parameters.at( 0 ) = _mvir0_;
		parameters.at( 1 ) = z();
		parameters.at( 2 ) = _c_;
		parameters.at( 3 ) = _tau_;
	}
	catch ( const std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}


/**
 *
 * @param parameter_names
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::get_parameter_names(std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	parameter_names.resize( num_parameters() );

	try
	{
		parameter_names.at( 0 ) = "mvir0";
		parameter_names.at( 1 ) = "z";
		parameter_names.at( 2 ) = "c";
		parameter_names.at( 3 ) = "tau";
	}
	catch ( const std::exception & )
	{
		return errorNOS( silent );
	}
	return 0;
}


/**
 *
 * @param f
 * @param silent
 * @return
 */
const int SALTSA::tNFW_profile::truncate_to_fraction( const double f,
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
		double new_tau_val = _taufm( f, bound(SMALL_FACTOR,std::fabs(0.1*(1-f)),0.00001) );
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

// SALTSA::point_mass_profile class methods
#if (1)

#if (1) // Constructors

/**
 *
 */
SALTSA::point_mass_profile::point_mass_profile()
{
	_mass_ = 0;
}

/**
 *
 * @param init_mass
 * @param init_z
 */
SALTSA::point_mass_profile::point_mass_profile( const double init_mass,
		const double init_z )
{
	set_mvir( init_mass );
	set_z( init_z );
}

#endif // End constructors

// Destructor
/**
 *
 */
SALTSA::point_mass_profile::~point_mass_profile()
{
}

#if (1) // Set functions

/**
 *
 * @param new_halo_mass
 * @param silent
 * @return
 */
const int SALTSA::point_mass_profile::set_mvir(
		const double new_halo_mass, bool silent )
{
	_mass_ = new_halo_mass;
	return 0;
}

/**
 *
 * @param num_parameters
 * @param parameters
 * @param silent
 * @return
 */
const int SALTSA::point_mass_profile::set_parameters(
		const unsigned int num_parameters,
		const std::vector< double > &parameters, bool silent )
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

/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::mvir() const
{
	return _mass_;
}
/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::mass() const
{
	return _mass_;
}

/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::mtot() const
{
	return _mass_;
}

/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::vvir() const
{
	return std::pow( 10 * Gc * H() * mvir(), 1. / 3. );
}
/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::rvir() const
{
	return vvir() / H() / 10;
}
/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::rs() const
{
	return 0;
}
/**
 *
 * @param silent
 * @return
 */
const double SALTSA::point_mass_profile::rt(
		const bool silent) const
{
	return 0;
}

/**
 *
 * @return
 */
const double SALTSA::point_mass_profile::hmvir() const
{
	return 0;
}

/**
 *
 * @param r
 * @return
 */
const double SALTSA::point_mass_profile::dens(
		const double r ) const
{
	double result = 0;

	result = ( r == 0 ? DBL_MAX : 0 );

	return result;
}
/**
 *
 * @param r
 * @param silent
 * @return
 */
const double SALTSA::point_mass_profile::enc_dens(
		const double r,
		const bool silent ) const
{
	return enc_mass( r ) / ( 4. / 3. * pi * cube( r ) );
}
/**
 *
 * @param r
 * @param silent
 * @return
 */
const double SALTSA::point_mass_profile::enc_mass(
		const double r,
		const bool silent ) const
{
	return _mass_;
}

/**
 *
 * @param parameters
 * @param silent
 * @return
 */
const int SALTSA::point_mass_profile::get_parameters( std::vector< double > & parameters,
		const bool silent ) const
{
	parameters.resize( num_parameters() );

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

/**
 *
 * @param parameter_names
 * @param silent
 * @return
 */
const int SALTSA::point_mass_profile::get_parameter_names( std::vector< std::string > & parameter_names,
		const bool silent ) const
{
	parameter_names.resize( num_parameters() );

	try
	{
		parameter_names.at( 0 ) = "mass";
		parameter_names.at( 1 ) = "z";
	}
	catch ( const std::exception & )
	{
		return errorNOS();
	}
	return 0;
}

/**
 *
 * @param f
 * @param silent
 * @return
 */
const int SALTSA::point_mass_profile::truncate_to_fraction( const double f,
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

// SALTSA::tfa_cache function implementations
/**
 *
 * @param in_params
 * @param out_params
 * @return
 */
const int SALTSA::tfa_cache::_calculate( const double in_params, double & out_params ) const
{
	try
	{
		out_params = -SALTSA::integrate_ltd( 0, SALTSA::zfa( in_params ) ) / c;
	}
	catch(const std::exception &e)
	{
		std::cerr << "ERROR: Could not calculate cache for " << _name_base() << "\n"
				<< "Exception: " << e.what() << "\n";
		std::cerr.flush();
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

// SALTSA::accel_function class methods
#if (1)

/**
 *
 * @param new_host
 * @return
 */
const int SALTSA::accel_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

/**
 *
 * @param in_param
 * @param out_param
 * @param silent
 * @return
 */
const int SALTSA::accel_function::operator()( const double & in_param,
double & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to accel_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = _host_ptr_->accel( in_param, silent );
	return 0;
}

/**
 *
 */
SALTSA::accel_function::accel_function()
{
	_host_ptr_ = NULL;
}
/**
 *
 * @param init_host
 */
SALTSA::accel_function::accel_function( const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end SALTSA::accel_function function implementations

// SALTSA::spherical_density_function class methods
#if (1)

/**
 *
 * @param new_host
 * @return
 */
const int SALTSA::spherical_density_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

/**
 *
 * @param in_param
 * @param out_param
 * @param silent
 * @return
 */
const int SALTSA::spherical_density_function::operator()(
		const double & in_param,
		double & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Host must be assigned to spherical_density_function before function can be called.\n";
		return NOT_SET_UP_ERROR;
	}
	out_param = 4 * pi * square( in_param )
			* _host_ptr_->dens( in_param );
	return 0;
}

/**
 *
 */
SALTSA::spherical_density_function::spherical_density_function()
{
	_host_ptr_ = NULL;
}
/**
 *
 * @param init_host
 */
SALTSA::spherical_density_function::spherical_density_function(
		const density_profile *init_host )
{
	set_host_ptr( init_host );
}

#endif // end SALTSA::spherical_density_function class methods

// SALTSA::solve_rhm_function class methods
#if (1)
/**
 *
 * @param new_host
 * @return
 */
const int SALTSA::solve_rhm_function::set_host_ptr(
		const density_profile *new_host )
{
	_host_ptr_ = new_host;
	return 0;
}

/**
 *
 * @param new_target_mass
 * @return
 */
const int SALTSA::solve_rhm_function::set_target_mass(
		const double new_target_mass )
{
	_target_mass_ = new_target_mass;
	return 0;
}

/**
 *
 * @param in_param
 * @param out_param
 * @param silent
 * @return
 */
const int SALTSA::solve_rhm_function::operator()( const double & in_param,
double & out_param, const bool silent ) const
{
	if ( _host_ptr_ == NULL )
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

/**
 *
 */
SALTSA::solve_rhm_function::solve_rhm_function()
{
	_host_ptr_ = NULL;
	_target_mass_ = 0;
}

/**
 *
 * @param init_host
 * @param new_target_mass
 */
SALTSA::solve_rhm_function::solve_rhm_function(
		const density_profile *init_host, const double new_target_mass )
{
	set_host_ptr( init_host );
	set_target_mass( new_target_mass );
}
#endif // end SALTSA::solve_rhm_function function implementations

#endif// end Class methods

/** Global Function Definitions **/
#if (1)

/**
 *
 * @param test_z
 * @return
 */
const double SALTSA::H( const double test_z )
{
	// Friedmann equation, assuming omega = -1
	if(test_z==0) return H_0;
	double zp1 = 1.+test_z;
	return H_0
			* std::sqrt( Omega_r * quart( zp1 )
							+ Omega_m * cube( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

/**
 *
 * @param a
 * @return
 */
const double SALTSA::zfa( const double a )
{
	return 1. / safe_d( a ) - 1.;
}
/**
 *
 * @param z
 * @return
 */
const double SALTSA::afz( const double z )
{
	return 1. / safe_d( 1 + z );
}

/**
 *
 * @param z
 * @return
 */
const double SALTSA::tfz( const double z )
{
	return SALTSA::tfa( afz( z ) );
}
/**
 *
 * @param a
 * @return
 */
const double SALTSA::tfa( const double a )
{
	return tfa_cache().get( a );
}
/**
 *
 * @param t
 * @return
 */
const double SALTSA::zft( const double t )
{
	return SALTSA::zfa( SALTSA::aft( t ) );
}
/**
 *
 * @param t
 * @return
 */
const double SALTSA::aft( const double t )
{
	SALTSA::tfa_cache cache;
	return cache.inverse_get( t );
}

// Functions to integrate out distances
#if (1)
/**
 *
 * @param z1
 * @param z2
 * @return
 */
const double SALTSA::integrate_ltd( const double z1, const double z2 )
{
	return SALTSA::integrate_distance( z1, z2, 3, 10000000 );
}
/**
 *
 * @param z
 * @return
 */
const double SALTSA::integrate_ltd( const double z )
{
	return SALTSA::integrate_distance( 0, z, 3, 10000000 );
}
/**
 *
 * @param z1_init
 * @param z2_init
 * @param mode
 * @param n
 * @return
 */
const double SALTSA::integrate_distance( const double z1_init,
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

	if(z1==z2) return 0;

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
				* std::sqrt( OR0 * inv_quart( a1 ) + OM0 * inv_cube( a1 )
								+ OK0 * inv_square( a1 ) + OL0 );
		OM = OM0 * square( H_0 / h1 ) * inv_cube( a1 );
		OR = OR0 * square( H_0 / h1 ) * inv_quart( a1 );
		OL = OL0 * square( H_0 / h1 );
		OK = 1 - OM - OR - OL;
	}

	HD = c / h1;

	for ( i = n; i >= 1; i-- )        // This loop is the numerical integration
	{
		a = ( i - 0.5 ) / n;              // Steadily decrease the scale factor
		// Comoving formula (See section 4 of Hogg, but I've added a radiation term too):
		adot = a * std::sqrt( OM * inv_cube(a) + OK * inv_square(a)
								+ OL + OR * inv_quart(a) ); // Note that "a" is equivalent to 1/(1+z)
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
		DM = ( 1 / std::sqrt( OK ) ) * sinh( std::sqrt( OK ) * DC );
	else if ( OK < -0.0001 )
		DM = ( 1 / std::sqrt( fabs( OK ) ) ) * sin( std::sqrt( fabs( OK ) ) * DC );
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

/**
 *
 * @param host
 * @param r
 * @param vr
 * @param vt
 * @return
 */
const double SALTSA::period( const SALTSA::density_profile *host,
		const double r, const double vr, const double vt )
{
	double mu = host->enc_mass( r ) * Gc;
	double v = quad_add( vr, vt );
	double a = -mu / 2 / safe_d( v * v / 2 - mu / safe_d( r ) );
	double result = (
			a > 0 ? 2 * pi * std::sqrt( cube(a) / mu ) : double( 0 ) );
	return result;
}

#endif // end other functions

#endif // end Global function definitions
