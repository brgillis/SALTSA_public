/**       @file SALTSA_example.cpp
 *
 *     Project: SALTSA_example
 *    Location: /SALTSA_example/SALTSA_example.cpp
 *
 *  Created on: 17 July 2014
 *      Author: brg
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <utility> // Just needed for some file loading here

#include "SALTSA.h"
#include "SALTSA_calculus.hpp" // This is just used to generate an orbit in this file.
#include "SALTSA_file_functions.h" // This is just used for file loading done here
#include "SALTSA_tSIS_profile.hpp" // In case we want to use this "user-defined" profile instead

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main( const int argc, const char *argv[] )
{

	// Set-up const values
#if (1) // This dummy compiler directive can be used for folding in Eclipse, and possibly other IDEs
	const int orbit_resolution = 1000000; // How many steps we take to integrate each orbit's path
	const int stripping_resolution = 1000; // Base number of steps to take to integrate stripping
	const int spline_points = 1000; // Number of points in the orbit we tell the stripping integrator
	                                // (It will use spline interpolation to estimate the rest.)
	const int spline_skip = SALTSA::max(orbit_resolution/spline_points,1); // How many points we skip between
	                                                                       // orbit points we report.

	// Toggle this flag to use a different mass profile, which is defined in header
	// SALTSA_tSIS_profile.hpp. This gives an example of how to use a user-defined profile with
	// SALTSA.
	const bool use_tSIS_profile=false;

	// Some information about the host and satellite haloes. Here we set this up to match the
	// set-up of Taylor and Babul 2004.

	const double host_z = 0, satellite_z = host_z; // Redshift
	const double host_c = 17.3, satellite_c = 17.3; // Concentration for truncated NFW profile

	// Mass of host and satellite. Note unit conversion here. SALTSA uses SI units at all times, so
	// to input something in units of Msun, we multiply it by unitconv::Msuntokg to get the value
	// in kg.
	const double host_mass = std::pow(10,14)*unitconv::Msuntokg;
	const double satellite_mass = std::pow(10,11)*unitconv::Msuntokg;

	// Set up density profiles representing the host now

	// Get a truncated NFW profile for the host, and a pointer to it (which we need to use to
	// exploit polymorphism).
	SALTSA::tNFW_profile host_group_tNFW(host_mass,host_z,host_c);

	// Or instead use a truncated SIS profile for the host
	SALTSA::tSIS_profile host_group_tSIS;

	host_group_tSIS.set_mvir(host_mass);
	host_group_tSIS.set_z(host_z);

	const SALTSA::density_profile *host_group;

	if(use_tSIS_profile)
		host_group = &host_group_tSIS;
	else
		host_group = &host_group_tNFW;

	const int num_host_parameters=host_group->num_parameters(); // How many parameters are used to define
	                                                            // the host halo.

	// Also get a point mass profile representing the group. We need this to properly calculate
	// the orbits, as Taylor and Babul 2004 used elliptical orbits, which we can only get with a
	// Keplerian potential (from a point mass).
	const SALTSA::point_mass_profile host_group_pm_val( host_mass, host_z );
	const SALTSA::density_profile *host_group_orbit = &host_group_pm_val; // And a pointer to it as well, again
	                                                                      // needed for polymorphism.

	// Get a truncated NFW profile for the satellite, and a pointer to it (which we need to use to
	// exploit polymorphism).
	SALTSA::tNFW_profile init_satellite_tNFW(satellite_mass,satellite_z,satellite_c);

	// Or instead use a truncated SIS profile for the satellite
	SALTSA::tSIS_profile init_satellite_tSIS;

	init_satellite_tSIS.set_mvir(host_mass);
	init_satellite_tSIS.set_z(host_z);

	const SALTSA::density_profile *init_satellite;

	if(use_tSIS_profile)
		init_satellite = &init_satellite_tSIS;
	else
		init_satellite = &init_satellite_tNFW;

	const int num_satellite_parameters=init_satellite->num_parameters(); // How many parameters are used to define
	                                                                     // the initial satellite.

	// The number of orbits we want the satellite to take
	const int num_periods = 3;

	// Get the total time, using the otvir method (orbital time at virial velocity), and the
	// step length for our orbit integration.
	const double total_time = num_periods*host_group->otvir(),
			t_step = total_time/orbit_resolution;

	// Set-up for the names of files we want to output
	const std::string output_file_name_base = "SALTSA_example_orbit_";
	const std::string output_file_name_tail = ".dat";

	// Set-up for orbit trajectories

	const double orbital_v_factor = 1; // How fast the satellite is at the virial radius compared to the
	                                   // virial velocity. 1 here means it starts at v = v_vir

	const int num_orbital_circs = 13; // Number of different circularities we want to test
	const double orbital_circularity[num_orbital_circs] = {1.00, .99, .95, .9, .8, .7, .6, .5, .4, .3, .2, .1, .05};
		// Array of the different circularities we'll be testing.

#endif // End set-up of const values

	// Declarations
#if (1)

	// Declarations for position, velocity, and time variables
	double x=0, y=0, z=0;
	double vx=0, vy=0, vz=0;
	double t=0;
	std::vector< double > host_parameters(num_host_parameters,0);

	SALTSA::accel_function accel_func(host_group_orbit);

	std::ofstream out;
	std::stringstream ss;

	std::string output_file_name = "";

	SALTSA::stripping_orbit test_orbit(host_group,init_satellite);

	// Vectors for setting up which parameters of the satellite we'll output
	std::vector<bool> satellite_output_parameters(num_satellite_parameters,false);
	std::vector<double> satellite_output_parameter_unitconvs(num_satellite_parameters,1);

	// Vectors for setting up which parameters of the host we'll output
	std::vector<bool> host_output_parameters(num_host_parameters,false);
	std::vector<double> host_output_parameter_unitconvs(num_host_parameters,1);

#endif // End declarations

	// Set-up for how to run SALTSA
#if (1)

	// Set up tuning parameters for SALTSA. Since we'll be doing a lot of orbits here, we override the
	// default tuning parameters. The "true" here causes these functions to also override the current
	// value of the tuning parameters when changing the default.
	test_orbit.set_default_tidal_stripping_amplification(0.625,true);
	test_orbit.set_default_tidal_stripping_deceleration(0.15,true);

	if(use_tSIS_profile)
	{
		// We only want to have sigma_v of the satellite be output here, so we set its value in the
		// output parameters array to true, and leave the rest as false
		satellite_output_parameters.at(0) = true; // Will output sigma_v

		// For the host, let's output all of its parameters
		host_output_parameters.at(0) = true; // Will output sigma_v
		host_output_parameters.at(1) = true; // Will output z

		// If we left it alone, we'd get the host sigma_v in m/s. But let's change that to km/s
		host_output_parameter_unitconvs.at(0) = unitconv::kmpstomps;
	}
	else
	{
		// We only want to have tau of the satellite be output here, so we set its value in the
		// output parameters array to true, and leave the rest as false
		satellite_output_parameters.at(3) = true; // Will output tau (and nothing else for the satellite)

		// For the host, let's output all of its parameters
		host_output_parameters.at(0) = true; // Will output mvir0
		host_output_parameters.at(1) = true; // Will output z
		host_output_parameters.at(2) = true; // Will output c
		host_output_parameters.at(3) = true; // Will output tau

		// If we left it alone, we'd get the host mass in kg. But let's change that to 10^10 Msun
		host_output_parameter_unitconvs.at(0) = unitconv::ttMsuntokg;
	}



	// Get the initial parameters array of the host. I'll use this here to demonstrate how you
	// would tell SALTSA about changes to the host over time.
	host_group->get_parameters(host_parameters);

#endif // set-up for how to run SALTSA

	// Orbit integration and stripping calculation
#if (1)

	// For each circularity
	for(int i = 0; i < num_orbital_circs; i += 1)
	{

        if(i!=0)
        {
        	// If this isn't the first orbit, we'll want to clear results from last time
        	test_orbit.clear();

        	// Set up the host and initial satellite pointers now
        	test_orbit = SALTSA::stripping_orbit(host_group, init_satellite);

        	// Alternatively, you could do the below:
        	// test_orbit.set_init_host(host_group);
        	// test_orbit.set_init_satellite(init_satellite);

        	// Or alternatively alternatively, if you want to use tNFW profiles, there's a built-in way
        	// to do that, without declaring an external profile:
        	// test_orbit.set_tNFW_init_host(host_mass,host_z,host_c);
        	// test_orbit.set_tNFW_init_satellite(satellite_mass,satellite_z,satellite_c);
        }

        // Let's tell the orbit what output parameters we want and what units we want them in.
        // Note that we don't tell it anything about the satellite's output parameter's units,
        // so they'll come out in default units. (In this case, only tau is output, which is unitless,
        // so that's what we want anyway.)
		test_orbit.set_satellite_output_parameters(num_satellite_parameters, satellite_output_parameters);
		test_orbit.set_host_output_parameters(num_host_parameters, host_output_parameters);
		test_orbit.set_host_parameter_unitconvs(num_host_parameters, host_output_parameter_unitconvs);

		// Set up the initial position of the orbit, starting at the virial radius
		x = 0;
		y = host_group_orbit->rvir();
		z = 0;

		// Determine the initial angle we want the satellite's velocity to be at relative to
		// a circular orbit
        double a = std::acos(orbital_circularity[i]/SALTSA::safe_d(orbital_v_factor));

		// Set up the initial velocity, depending on the angle calculated above
		vx = host_group_orbit->vvir()*orbital_v_factor*cos(a);
		vy = -host_group_orbit->vvir()*orbital_v_factor*sin(a);
		vz = 0;

		// Start at time zero (in seconds)
		t = 0;

		// Set up the name of the file we want to output to
		ss.str("");
		ss << output_file_name_base << i << output_file_name_tail;
		output_file_name = ss.str();

		// Tell the user we've started to work on this orbit
		std::cout << "Working on " << output_file_name << "... " << std::flush;

		// Now we'll integrate out the orbit. We'll use a leapfrog(-ish) integrator I programmed
		// for this.
		for( int j=0; j < orbit_resolution; j++)
		{

			// Every spline_skip point, we'll tell the stripping_orbit about this position
			if( SALTSA::divisible(j,spline_skip) )
			{
				// Tell SALTSA about this point. We won't always get these values from orbit integration -
				// they could also be extracted from a file, for instance.
				test_orbit.add_point(x, y, z, vx, vy, vz, t);

				// What if we didn't know v? In that case, we could do:
				// test_orbit.add_point(x, y, z, t);
				// SALTSA would then differentiate the position over time and apply a smoothing
				// kernel to estimate v. It's not perfect though - we have to smooth to keep from
				// amplifying noise, but this tends to decrease the overall velocity variation. This
				// can be tuned a bit if it has to be used to balance between getting rid of noise
				// and not suppressing velocity variation too much.

				// Let's also tell the stripping_orbit about the state of the host at this time.
				// Here, we aren't actually varying it, but you would do the same thing if it were
				// varying.

				test_orbit.add_host_parameter_point(host_parameters, t);
			}

			// Integrate to next point with leapfrog-ish method (Not strictly leapfrog since we don't
			// start out with x and v staggered, but it converges to the same result for large N).
			SALTSA::leapfrog_step( x, y, z, vx, vy, vz, t_step, &accel_func );
			t += t_step;
		}

		// What if the satellite changed hosts at some point? Then you'd add in all the points as
		// they were for the current host, and do the following command to tell SALTSA to split
		// the orbit up into discrete segments: (try uncommenting it)
		// test_orbit.add_discontinuity_time(total_time/3.);

		// Now, let's calculate and output the data for this orbit

		test_orbit.set_record_full_data(true); // Tell SALTSA to record detailed output data when it's run
		                                       // If this isn't done, it may end up running a second time
		                                       // when it finds out it needs detailed data.
		test_orbit.set_resolution(stripping_resolution); // Tell SALTSA the resolution to use for
		                                                 // calculating stripping

		// Tell SALTSA to calculate now. This doesn't have to be called explicitly - it will be done
		// when needed and the results cached for quick access afterward.
		if(test_orbit.calc())
		{
			std::cerr << "Something went wrong in calculating stripping.\n";
			return 1;
		}

		// Output to the appropriate file
		out.open(output_file_name.c_str()); // Open the file
		test_orbit.print_full_data( &out ); // Tell SALTSA to output to it
		out.close(); // Close the file
		out.clear(); // Clear the stream (needed on certain systems when it'll be reused)

		std::cout << "Done!\n";

		// And that's it for this orbit. We'll go through all of them and calculate stripping for
		// various circularities, storing the results in ASCII text files.
	}
#endif

	// File reading and stripping calculation
#if (1)
	// Now, let's demonstrate how MATSA can compare to another stripping estimate.
	// We'll start by loading the data we generated before, and then doing various comparison
	// runs, so we can see how the results change with different integration types.

	// Set up
#if (1)
	// Get the name of the file we'll use

	// We'll use the orbit with circularity .8 for comparison
	unsigned int comparison_index = 4;

	ss.str("");
	ss << output_file_name_base << comparison_index << output_file_name_tail;
	std::string comparison_file_name = ss.str();

	// Declare vectors to store the orbit data in
	std::vector<double> t_data;

	std::vector<double> x_data;
	std::vector<double> y_data;
	std::vector<double> z_data;

	std::vector<double> vx_data;
	std::vector<double> vy_data;
	std::vector<double> vz_data;

	std::vector<double> host_mvir0_data;
	std::vector<double> host_sigma_v_data;
	std::vector<double> host_z_data;
	std::vector<double> host_c_data;
	std::vector<double> host_tau_data;

	std::vector<double> comparison_mret_data;
#endif // Set up

	// Loading the data
#if (1)
	// Since the data has been saved in a standard ASCII table, we can use the loading
	// functions I've programmed to load it relatively easily. Since the saved data
	// can have a variable number of columns, we'll use a function which parses the
	// header to find the right column for each vector. (Try varying which satellite
	// parameters are output to test this. The number of those output will affect the
	// indices of the host parameters, but this function will pick them out appropriately
	// anyway. This will only fail if the density profile for the host is changed, which
	// would change the names of its parameters.)

	// To use this function, we'll need to set up "key"s, which connect strings to look
	// for in the header to vectors to put data into. We'll then pass a vector of these
	// keys to the load_table_columns function.

	// Declare the key vector
	typedef std::pair<std::string, std::vector<double> * > key;
	std::vector<key> header_keys(0);

	// Load it with keys for each column we want to use, using std::make_pair to easily
	// construct keys
	header_keys.push_back( std::make_pair("t",&t_data) );

	header_keys.push_back( std::make_pair("x",&x_data) );
	header_keys.push_back( std::make_pair("y",&y_data) );
	header_keys.push_back( std::make_pair("z",&z_data) );

	header_keys.push_back( std::make_pair("vx",&vx_data) );
	header_keys.push_back( std::make_pair("vy",&vy_data) );
	header_keys.push_back( std::make_pair("vz",&vz_data) );

	if(use_tSIS_profile)
	{
		header_keys.push_back( std::make_pair("Host_sigma_v",&host_sigma_v_data) );
		header_keys.push_back( std::make_pair("Host_z",&host_z_data) );

		header_keys.push_back( std::make_pair("m_ret",&comparison_mret_data) );
	}
	else
	{
		header_keys.push_back( std::make_pair("Host_mvir0",&host_mvir0_data) );
		header_keys.push_back( std::make_pair("Host_z",&host_z_data) );
		header_keys.push_back( std::make_pair("Host_c",&host_c_data) );
		header_keys.push_back( std::make_pair("Host_tau",&host_tau_data) );

		header_keys.push_back( std::make_pair("m_ret",&comparison_mret_data) );
	}

	// And now actually load it
	// An exception will be thrown here if a column we want isn't found in the file, or if it's
	// improperly formatted.
	SALTSA::load_table_columns(comparison_file_name, header_keys);

	// Now, all the vectors should be loaded. Let's pass the data into SALTSA

	test_orbit.clear();
	test_orbit = SALTSA::stripping_orbit(host_group, init_satellite);

	// We'll simulate the original data through skipping some points (since the data we read
	// reflects the adaptive step size)
	double t_skip_interval = total_time/spline_points/2; // /2 correction
	double t_next = 0;

	for(unsigned int i=0; i<t_data.size(); i++)
	{
		// Load the first point after each t interval (to simulate the way we initially
		// assigned points)
		if(t_data.at(i)*unitconv::Gyrtos >= t_next)
		{
			test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
					y_data.at(i)*unitconv::kpctom,
					z_data.at(i)*unitconv::kpctom,
					vx_data.at(i)*unitconv::kmpstomps,
					vy_data.at(i)*unitconv::kmpstomps,
					vz_data.at(i)*unitconv::kmpstomps,
					t_data.at(i)*unitconv::Gyrtos,
					comparison_mret_data.at(i));
			host_parameters.clear();
			if(use_tSIS_profile)
			{
				host_parameters.push_back(host_sigma_v_data.at(i)*unitconv::kmpstomps);
				host_parameters.push_back(host_z_data.at(i));
			}
			else
			{
				host_parameters.push_back(host_mvir0_data.at(i)*unitconv::ttMsuntokg);
				host_parameters.push_back(host_z_data.at(i));
				host_parameters.push_back(host_c_data.at(i));
				host_parameters.push_back(host_tau_data.at(i));
			}

			test_orbit.add_host_parameter_point(host_parameters,
					t_data.at(i)*unitconv::Gyrtos);

			t_next =  t_data.at(i)*unitconv::Gyrtos+t_skip_interval;
		}
	}

#endif // Loading the data

	// Base comparison to see how far off this dataset is just through numeric errors
#if (1)

	// General set-up for the orbit
	test_orbit.set_satellite_output_parameters(num_satellite_parameters, satellite_output_parameters);
	test_orbit.set_host_output_parameters(num_host_parameters, host_output_parameters);
	test_orbit.set_host_parameter_unitconvs(num_host_parameters, host_output_parameter_unitconvs);
	test_orbit.set_resolution(stripping_resolution);

	test_orbit.set_resolution(stripping_resolution);

	test_orbit.set_record_full_data(true);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_base_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	double diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done!\n";
	std::cout << "Base difference is " << diff << ". Differences of this scale should be considered negligible.\n";

#endif // Base comparison

	// Testing changes to integration and set-up parameters
#if (1)

	// Let's try changing the resolution
#if (1)

	// High resolution
	test_orbit.set_resolution(10*stripping_resolution);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_hires_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will automatically calculate here
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Low resolution
	test_orbit.set_resolution(0.1*stripping_resolution);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_lores_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change resolution back to how it was before
	test_orbit.set_resolution(stripping_resolution);

#endif // Let's try changing the resolution

	// Different interpolation types
#if (1)

	// Linear interpolation
	test_orbit.set_interpolation_type(SALTSA::stripping_orbit::LINEAR);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_lininterp_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// No interpolation (LOWER just uses the closest value below it)
	test_orbit.set_interpolation_type(SALTSA::stripping_orbit::LOWER);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_nointerp_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change interpolation back to how it was before
	test_orbit.set_interpolation_type(SALTSA::stripping_orbit::UNSET); // Which results in spline interpolation here

#endif // Different interpolation types

	// Different adaptive step size methods
#if (1)

	// Constant-time step size
	test_orbit.set_step_length_power(0.);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_noadaptstep_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Strongly adaptive step size
	test_orbit.set_step_length_power(3.);
	test_orbit.set_step_factor_min(0.001);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_strongadaptstep_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Note that the strongly adaptive step size takes many more steps (about 80% as many as the
	// high-resolution run), but it is indeed more accurate. It comes to the same result as the
	// high-resolution run, but with a bit fewer steps. So it's likely more efficient to make
	// the step size more strongly adaptive than to increase the resolution directly, as this
	// results in more points near pericentre, where they're most needed.

	// Change adaptive step length back to how it was before
	test_orbit.reset_step_length_power();
	test_orbit.reset_step_factor_min();

#endif // Different adaptive step size methods

	// Skipping velocity values
#if (1)

	// For this one, we'll actually have to go back and completely reload the orbit to show
	// how this can be done

	// Let's only give velocity for every other point
	bool record_velocity_for_this_point = true;

	ss.str("");
	ss << output_file_name_base << comparison_index << "_skipsomev_comp" << output_file_name_tail;
	output_file_name = ss.str();

	test_orbit.clear_points(); // This clears all position, velocity, and mass comparison points,
	                           // but not host parameter points

	t_next = 0;

	for(unsigned int i=0; i<t_data.size(); i++)
	{
		if(t_data.at(i)*unitconv::Gyrtos >= t_next)
		{
			if(record_velocity_for_this_point)
			{
				test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
						y_data.at(i)*unitconv::kpctom,
						z_data.at(i)*unitconv::kpctom,
						vx_data.at(i)*unitconv::kmpstomps,
						vy_data.at(i)*unitconv::kmpstomps,
						vz_data.at(i)*unitconv::kmpstomps,
						t_data.at(i)*unitconv::Gyrtos,
						comparison_mret_data.at(i));
				record_velocity_for_this_point = false;
			}
			else
			{
				test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
						y_data.at(i)*unitconv::kpctom,
						z_data.at(i)*unitconv::kpctom,
						t_data.at(i)*unitconv::Gyrtos,
						comparison_mret_data.at(i));
				record_velocity_for_this_point = true;
			}
			t_next = t_data.at(i)*unitconv::Gyrtos + t_skip_interval;
		}
	}

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Let's not give any velocity points
	ss.str("");
	ss << output_file_name_base << comparison_index << "_skipallv_comp" << output_file_name_tail;
	output_file_name = ss.str();

	test_orbit.clear_points(); // This clears all position, velocity, and mass comparison points,
	                           // but not host parameter points

	t_next = 0;

	for(unsigned int i=0; i<t_data.size(); i++)
	{
		if(t_data.at(i)*unitconv::Gyrtos >= t_next)
		{
			test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
					y_data.at(i)*unitconv::kpctom,
					z_data.at(i)*unitconv::kpctom,
					t_data.at(i)*unitconv::Gyrtos,
					comparison_mret_data.at(i));
			t_next = t_data.at(i)*unitconv::Gyrtos + t_skip_interval;
		}
	}

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Not so bad with a pretty circular orbit, like we're using here, but try this also with a
	// more eccentric orbit (like number 9), and see how much worse it is in that case.

	// Lesson here: For particularly eccentric orbits, the differentiation performed by SALTSA
	// to estimate velocity tends to systematically underestimate the spikes near pericentre,
	// resulting in an underestimate of tidal stripping there compared to if velocity data
	// is given.

#endif // Skipping velocity values

	// Telling SALTSA about fewer points on the orbit
#if (1)
	// Let's see what happens if we only tell SALTSA about one in every ten points
	// on the orbit.

	// Once more, we'll have to reload the orbit

	// Let's only load one tenth as many points
	t_skip_interval *= 10;
	t_next = 0;

	ss.str("");
	ss << output_file_name_base << comparison_index << "_fewerorbitpoints_comp" << output_file_name_tail;
	output_file_name = ss.str();

	test_orbit.clear_points(); // This clears all position, velocity, and mass comparison points,
	                           // but not host parameter points
	t_next = 0;

	for(unsigned int i=0; i<t_data.size(); i++)
	{
		// Load the first point after each t interval (to simulate initially using
		// only one-tenth the points. If we just only input every tenth, we'll also
		// mimic the adaptive step size used, which we don't want to do here.)
		if(t_data.at(i)*unitconv::Gyrtos >= t_next)
		{
			test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
					y_data.at(i)*unitconv::kpctom,
					z_data.at(i)*unitconv::kpctom,
					vx_data.at(i)*unitconv::kmpstomps,
					vy_data.at(i)*unitconv::kmpstomps,
					vz_data.at(i)*unitconv::kmpstomps,
					t_data.at(i)*unitconv::Gyrtos,
					comparison_mret_data.at(i));
			t_next = t_data.at(i)*unitconv::Gyrtos + t_skip_interval;
		}
	}

	t_skip_interval /= 4;

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out );
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// (Some ringing occurs here in the comparison m_ret data; this is due to spline
	// interpolation being used to guess intermediate points so a value can be plotted
	// for them.)

	// This looks bad... can we compensate for this by increasing resolution?

	test_orbit.set_resolution(10*stripping_resolution);

	ss.str("");
	ss << output_file_name_base << comparison_index << "_fewerorbitpoints_hires_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will automatically calculate here
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	test_orbit.set_resolution(stripping_resolution);

	// Lesson here: For eccentric orbits, having too few actual data points along the orbit's
	// path can lead to an underestimation of stripping, which can't simply be fixed with
	// increasing resolution in SALTSA - only more orbit snapshots will help.
	//
	// This could be a significant issue, as we have to use 1000 snapshots here to get
	// decent results, but most simulations will only give on order of 100 snapshots
	// (where we see a problem).
	//
	// Future versions of SALTSA may try to account for this, for instance by using the
	// orbit's circularity to better estimate how close its pericentre is.

#endif // Telling SALTSA about fewer points on the orbit

#endif // Testing changes to integration and set-up parameters

	// Testing changes to tuning parameters
#if (1)

	// Let's first reset the orbit and load it up as before

	test_orbit.clear_points();
	t_next = 0;

	for(unsigned int i=0; i<t_data.size(); i++)
	{
		// Load the first point after each t interval (to simulate initially using
		// only one-tenth the points. If we just only input every tenth, we'll also
		// mimic the adaptive step size used, which we don't want to do here.)
		if(t_data.at(i)*unitconv::Gyrtos >= t_next)
		{
			test_orbit.add_point(x_data.at(i)*unitconv::kpctom,
					y_data.at(i)*unitconv::kpctom,
					z_data.at(i)*unitconv::kpctom,
					vx_data.at(i)*unitconv::kmpstomps,
					vy_data.at(i)*unitconv::kmpstomps,
					vz_data.at(i)*unitconv::kmpstomps,
					t_data.at(i)*unitconv::Gyrtos,
					comparison_mret_data.at(i));
			t_next = t_data.at(i)*unitconv::Gyrtos + t_skip_interval;
		}
	}

	// We'll look at each tuning parameter in turn

	// Tuning stripping amplification (A_s in the paper)
#if (1)

	// Increase the amplification by 0.2
	test_orbit.set_tidal_stripping_amplification(0.2 +
			test_orbit.default_tidal_stripping_amplification());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_histripamp_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change stripping amplification back to how it was before
	test_orbit.reset_tidal_stripping_amplification();

#endif // Tuning stripping amplification

	// Tuning stripping deceleration (alpha_s in the paper)
#if (1)

	// Increase the deceleration by 0.2
	test_orbit.set_tidal_stripping_deceleration(0.2 +
			test_orbit.default_tidal_stripping_deceleration());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_histripdecel_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change stripping deceleration back to how it was before
	test_orbit.reset_tidal_stripping_deceleration();

#endif // Tuning stripping deceleration

	// Tuning shocking amplification (A_h in the paper)
#if (1)

	// Increase the amplification by 0.5
	test_orbit.set_tidal_shocking_amplification(0.5 +
			test_orbit.default_tidal_shocking_amplification());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_hishockamp_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change shocking amplification back to how it was before
	test_orbit.reset_tidal_shocking_amplification();

#endif // Tuning shocking amplification (A_h in the paper)

	// Tuning shocking power (alpha_h in the paper)
#if (1)

	// Increase the power by 1
	test_orbit.set_tidal_shocking_power(1 +
			test_orbit.default_tidal_shocking_power());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_hishockpow_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Change shocking power back to how it was before
	test_orbit.reset_tidal_shocking_power();

#endif // Tuning shocking power

	// Tuning shocking persistance
#if (1)

	// Increase the persistance by 1
	test_orbit.set_tidal_shocking_persistance(1 +
			test_orbit.default_tidal_shocking_persistance());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_hishockpersist_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// This looks very similar to changing amplification. Can they cancel each other out?

	// Decrease the amplification by 0.5
	test_orbit.set_tidal_shocking_amplification(-0.5 +
			test_orbit.default_tidal_shocking_amplification());

	ss.str("");
	ss << output_file_name_base << comparison_index << "_hishockpersistloshockamp_comp" << output_file_name_tail;
	output_file_name = ss.str();

	out.open(output_file_name.c_str());

	std::cout << "Working on " << output_file_name << "... " << std::flush;
	test_orbit.print_full_data( &out ); // It will recalculate now, with lower resolution
	diff = test_orbit.quality_of_fit(); // Get a measurement for how different the two curves are
	std::cout << "Done! Difference measure: " << diff << "\n";

	out.close();
	out.clear();

	// Lesson here: Shocking persistance and amplification are pretty degenerate, so it's only
	// worth tuning one of them. They also require much larger changes in their tuning values to
	// have significant effects on the orbits than the stripping tuning values do. It's likely
	// worth implementing a prior on them in any fitting to avoid them being fit to extreme
	// values.
	//
	// We can do a similar comparison to show that shocking power is also degenerate with these.
	// I recommend only tuning shocking amplification.

	// Change shocking persistance and amplification back to how they were before
	test_orbit.reset_tidal_shocking_persistance();
	test_orbit.reset_tidal_shocking_amplification();



#endif // Tuning shocking power

#endif // Testing changes to tuning parameters

#endif // File reading and stripping calculation

	return 0;
}
