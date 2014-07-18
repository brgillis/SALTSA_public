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
#include <sstream>
#include <fstream>

#include "SALTSA.h"
#include "SALTSA_calculus.hpp" // This is just used to generate an orbit in this file.

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
	const int stripping_resolution = 100000; // Base number of steps to take to integrate stripping
	const int spline_points = 100000; // Number of points in the orbit we tell the stripping integrator
	                                // (It will use spline interpolation to estimate the rest.)
	const int spline_skip = SALTSA::max(orbit_resolution/spline_points,1); // How many points we skip between
	                                                                       // orbit points we report.

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
	const SALTSA::tNFW_profile host_group_val( host_mass, host_z, host_c),
			*host_group = &host_group_val;
	const int num_host_parameters=host_group->num_parameters(); // How many parameters are used to define
	                                                            // the host halo.

	// Also get a point mass profile representing the group. We need this to properly calculate
	// the orbits, as Taylor and Babul 2004 used elliptical orbits, which we can only get with a
	// Keplerian potential (from a point mass).
	const SALTSA::point_mass_profile host_group_pm_val( host_mass, host_z );
	const SALTSA::density_profile *host_group_orbit = &host_group_pm_val; // And a pointer to it as well, again
	                                                                      // needed for polymorphism.

	// And a density profile and pointer to it representing the initial satellite.
	const SALTSA::tNFW_profile init_satellite_val( satellite_mass, satellite_z, satellite_c),
			*init_satellite = &init_satellite_val;
	const int num_satellite_parameters=init_satellite->num_parameters(); // How many parameters are used to define
	                                                                     // the initial satellite.

	// The number of orbits we want the satellite to take
	const int num_periods = 3;

	// Get the total time, using the otvir method (orbital time at virial velocity), and the
	// step length for our orbit integration.
	const double total_time = num_periods*host_group->otvir(),
			t_step = total_time/orbit_resolution;

	// Set-up for the names of files we want to output
	const std::string orbital_file_name_base = "SALTSA_example_orbit_";
	const std::string orbital_file_name_tail = ".dat";

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

	SALTSA::phase current_phase(x,y,z,vx,vy,vz,t);
	SALTSA::accel_function accel_func(host_group_orbit);

	std::ofstream out;
	std::stringstream ss;

	std::string orbital_file_name = "";

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
	test_orbit.set_default_tidal_stripping_amplification(0.65,true);
	test_orbit.set_default_tidal_stripping_deceleration(0.125,true);

	// We only want to have tau of the satellite be output here, so we set its value in the
	// output parameters array to true, and leave the rest as false
	satellite_output_parameters.at(3) = true; // Will output tau (and nothing else for the satellite)

	// For the host, let's output all of its parameters
	host_output_parameters.at(0) = true; // Will output mvir0
	host_output_parameters.at(1) = true; // Will output z
	host_output_parameters.at(2) = true; // Will output c
	host_output_parameters.at(3) = true; // Will output tau

	// If we left it alone, we'd get the host mass in kg. But let's change that to 10^10 Msun
	host_output_parameter_unitconvs.at(0) = unitconv::ttMsuntokg; // Will output mvir0 in 10^10 Msun

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
		ss << orbital_file_name_base << i << orbital_file_name_tail;
		orbital_file_name = ss.str();

		// Tell the user we've started to work on this orbit
		std::cout << "Working on " << orbital_file_name << "... " << std::flush;

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
		out.open(orbital_file_name.c_str()); // Open the file
		test_orbit.print_full_data( &out ); // Tell SALTSA to output to it
		out.close(); // Close the file
		out.clear(); // Clear the stream (needed on certain systems when it'll be reused)

		std::cout << "Done!\n";

		// And that's it for this orbit. We'll go through all of them and calculate stripping for
		// various circularities, storing the results in ASCII text files.
	}
#endif

	return 0;
}
