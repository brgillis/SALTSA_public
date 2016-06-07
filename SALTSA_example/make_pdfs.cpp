#include <make_pdfs.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
#include "brg/physics/density_profile/tNFW_profile.h"
#include "brg/file_access/open_file.hpp"
#include "brg/physics/SALTSA/stripping_orbit.h"
#include "brg/math/random/random_functions.hpp"

using namespace std;
namespace unitconv = brgastro::unitconv;

int main( int argc, char* argv[] )
{

	/*******************************Adjustable Parameters******************************/
	const char* datafile = argv[1];

	const bool use_interlopers = true;
	const char* interloperfile = { '\0' };
	if ( use_interlopers )
		interloperfile = argv[3];

	const int nsteps = 42;
	const int Nbinsr = 100;
	const int Nbinsv = 100;

	const float Rmin = 0;
	const float Rmax = 2.5;
	const float Vmin = 0;
	const float Vmax = 2.0;
	const float mrmin = 0;
	const float mrmax = 1;

	const unsigned int orbit_chunk_size = 50000;

	const float mass_at_infall_cut = (float)1.2411378667E10; //=10^11.9/64; //=10^12.9/64 //M_sun physical (not h^-1)

	const bool use_infalling = false; //use halos which default to infall a=1.0011 by not finding another robust time

	const bool quiet = true; //explicit warning messages for dropped halos?

	const bool record_all_orbits = false; // Record detailed data for each orbit. Overrides record_random_orbits if set to true
	const bool record_random_orbits = true; // Record a random subsample of orbits
	const double record_random_chance = 0.0001; // 0.01% chance for a given orbit to be recorded

	const bool debug_info = false;

	const bool record_general_orbit_data = true; // Records basic info about each orbit in a table
	const string general_orbit_data_file_name = "general_orbit_data.dat";

	const string output_file_name_base = "orbit_data_";
	const string output_file_name_tail = ".dat";

	brgastro::stripping_orbit::set_default_tidal_stripping_amplification(0.825);
	brgastro::stripping_orbit::set_default_tidal_stripping_deceleration(0.025);
	brgastro::stripping_orbit::set_default_tidal_shocking_amplification(3.0);
	brgastro::stripping_orbit::set_default_tidal_shocking_power(-2.5);

	//const int max_num_threads = 1; // Number of parallel processors to use

	//omp_set_num_threads( max_num_threads ); // Limit to maximum threads

	/**********************************************************************************/

	//create PDF array and initialize all values to 0
	int*** pdfs = new int**[Nbinsv];
	for ( int oloop = 0; oloop < Nbinsv; oloop++ )
	{
		pdfs[oloop] = new int*[Nbinsr];
		for ( int loop = 0; loop < Nbinsr; loop++ )
		{
			pdfs[oloop][loop] = new int[nsteps];
			for ( int iloop = 0; iloop < nsteps; iloop++ )
				pdfs[oloop][loop][iloop] = 0;
		}
	}

	//define bins in R-V space

	float Rbins[Nbinsr + 1];
	float Vbins[Nbinsv + 1];

	for ( int loop = 0; loop <= Nbinsr; loop++ )
	{
		Rbins[loop] = Rmin + ( Rmax - Rmin ) * loop * 1.0 / Nbinsr; //avoid integer division
	}
	for ( int loop = 0; loop <= Nbinsv; loop++ )
	{
		Vbins[loop] = Vmin + ( Vmax - Vmin ) * loop * 1.0 / Nbinsv; //avoid integer division
	}

	//define scale factor bins

	float mrbins[nsteps + 1];
	float mrs[nsteps];
	for ( int i = 0; i < nsteps; i++ )
		mrs[i] = mrmin + ( mrmax - mrmin ) * i / ( nsteps - 1 );
	for ( int loop = 0; loop < nsteps + 1; loop++ )
	{
		if ( loop == 0 )
		{
			mrbins[loop] = -.5 * ( mrs[loop + 1] - mrs[loop] ) + mrs[loop];
		}
		else if ( loop == nsteps )
		{
			mrbins[loop] = .5 * ( mrs[loop - 1] - mrs[loop - 2] )
					+ mrs[loop - 1];
		}
		else
		{
			mrbins[loop] = -.5 * ( mrs[loop] - mrs[loop - 1] ) + mrs[loop];
		}
	}

	ofstream outfile;

	const char* ofilename = argv[2];

	outfile.open( ofilename, ofstream::out | ofstream::trunc );

	outfile << Nbinsr << endl;
	outfile << Nbinsv << endl;
	outfile << nsteps << endl;
	for ( int loop = 0; loop < Nbinsr + 1; loop++ )
		outfile << Rbins[loop] << " ";
	outfile << endl;
	for ( int loop = 0; loop < Nbinsv + 1; loop++ )
		outfile << Vbins[loop] << " ";
	outfile << endl;
	for ( int loop = 0; loop < nsteps + 1; loop++ )
		outfile << mrbins[loop] << " ";
	outfile << endl;

	int dropped = 0;
	int interloper_count = 0;
	int output_file_count = 0;

	// Set up output file for general orbit data
	ofstream general_orbit_data_file;
	if(record_general_orbit_data)
	{
		// Try to open the file
		brgastro::open_file(general_orbit_data_file,general_orbit_data_file_name);

		// Output the header
		general_orbit_data_file << "# ID\tfirst_infall_mass_ttMsun\tfinal_host_mass_ttMsun\t"
				<< "r_kpc\tr/rvir\tR_kpc\tR/rvir\t"
				<< "v_km_per_s\tv/vvir\tv_los_km_per_s\tv_los/vvir\t"
				<< "t_first_infall_Gyr\tt_last_infall_Gyr\t"
				<< "z_first_infall\tz_last_infall\t"
				<< "frac_M_retained\n";
	}

	if ( use_interlopers == true )
	{
		vector< Interloper > Interlopers;
		cout << "Reading interloper data... ";
		cout.flush();
		read_interloper_data( interloperfile, Interlopers );
		cout << "Read in " << Interlopers.size() << " interlopers." << endl;
		int drop_bins = 0;
		int drop_mass = 0;

		for ( unsigned int loop = 0; loop < Interlopers.size(); loop++ )
		{
			//locate the correct bin for this interloper
			bool found[2] = { false, false };
			int r;
			int v;
			int a = nsteps - 1;
			//mass check
			if ( Interlopers[loop].get_mass() < mass_at_infall_cut )
			{
				drop_mass++;
				if ( quiet == false )
					cout << "Cut an interloper for low mass." << endl;
				continue;
			}
			// Output it to the general orbit data file
			if(record_general_orbit_data)
			{
				brgastro::tNFW_profile host(Interlopers[loop].get_host_mass()*unitconv::Msuntokg,0);
				general_orbit_data_file << loop << "\t" << Interlopers[loop].get_mass()*unitconv::Msuntokg*unitconv::kgtottMsun <<
						"\t" << host.mvir()*unitconv::kgtottMsun << "\t" <<
						-1 << "\t" <<
						-1 << "\t" <<
						Interlopers[loop].get_R()*host.rvir()*unitconv::mtokpc <<
						"\t" << Interlopers[loop].get_R() << "\t" <<
						-1 << "\t" <<
						-1 << "\t" <<
						fabs(Interlopers[loop].get_V())*host.vvir()*unitconv::mpstokmps <<
						"\t" << fabs(Interlopers[loop].get_V()) << "\t" <<
						0 << "\t" <<
						0 << "\t" <<
						0 << "\t" <<
						0 << "\t" <<
						1 << "\n";
			}
			//find V bin
			for ( int vloop = 0; vloop < Nbinsv; vloop++ )
			{
				if ( Interlopers[loop].get_V() > Vbins[vloop]
						&& Interlopers[loop].get_V() <= Vbins[vloop + 1] )
				{
					v = vloop;
					found[0] = true;
					break;
				}
			}
			v = Nbinsv-1;
			//find R bin
			for ( int rloop = 0; rloop < Nbinsr; rloop++ )
			{
				if ( Interlopers[loop].get_R() > Rbins[rloop]
						&& Interlopers[loop].get_R() <= Rbins[rloop + 1] )
				{
					found[1] = true;
					break;
				}
				r = rloop;
			}
			r = Nbinsr-1;
			if ( found[0] == true && found[1] == true )
			{
				//increment interloper bin
				interloper_count++;
				pdfs[v][r][a]++;
				continue;
			}
			else
			{
				drop_bins++;
				if ( quiet == false )
					cout << "No R-V bin found for interloper R="
							<< Interlopers[loop].get_R() << " V="
							<< Interlopers[loop].get_V() << "." << endl;
				continue;
			}
		}
		cout << "Dropped " << drop_bins
				<< " interlopers for no bin found, and " << drop_mass
				<< " interlopers for mass cut." << endl;
		cout << "Total interlopers used: " << interloper_count << "." << endl;
	}

	vector< Orbit > Orbits;

	unsigned int orbit_N_start = 0;

	do {
		cout << "Reading orbit data... ";
		cout.flush();
		read_orbit_data( datafile, Orbits, orbit_N_start, orbit_chunk_size );
		orbit_N_start += orbit_chunk_size;
		cout << "Read in " << Orbits.size() << " orbits." << endl;

		//loop over each orbit
		#pragma omp parallel for
		for ( unsigned int loop = 0; loop < Orbits.size(); loop++ )
		{
			if(debug_info)
			{
				std::cout << "Looking at orbit " << loop << endl;
			}

			Orbit current_orbit;

			#pragma omp critical
			{
				current_orbit = Orbits[loop];
			}

			//locate the correct bin for this orbit
			int r = Nbinsr - 1;
			int v = Nbinsv - 1;
			int mr = nsteps - 1;
			bool found[4] = { false, false, false, false };
			for ( int rloop = 0; rloop < Nbinsr; rloop++ )
			{
				int t_now = current_orbit.get_nsteps() - 1;
				if ( current_orbit.rfproj( t_now, 0 ) > Rbins[rloop]
						&& current_orbit.rfproj( t_now, 0 ) <= Rbins[rloop + 1] )
				{
					r = rloop;
					found[0] = true;
					break;
				}
			}
			for ( int vloop = 0; vloop < Nbinsv; vloop++ )
			{
				int t_now = current_orbit.get_nsteps() - 1;
				if ( abs( current_orbit.vfproj( t_now, 0 ) ) > Vbins[vloop]
						&& abs( current_orbit.vfproj( t_now, 0 ) )
								<= Vbins[vloop + 1] )
				{
					v = vloop;
					found[1] = true;
					break;
				}
			}

			if ( found[0] == true && found[1] == true )
			{
				double fmret = -1; // Flag to show it hasn't been calculated yet
				//compute infall time for this orbit
				int t_now = current_orbit.get_nsteps() - 1;
				int mrloop_at_infall = 0;
				for ( int mrloop = t_now - 1; mrloop >= 0; mrloop-- )
				{
					if ( current_orbit.get_parent( mrloop ) == 1
							&& current_orbit.get_fppb( mrloop ) == 1 )
					{
						if ( current_orbit.get_snp_mass( mrloop )
								== current_orbit.get_fppb_mass( mrloop ) )
						{
							if ( ( current_orbit.get_snp_mass( mrloop ) > 0 )
									&& ( current_orbit.get_fppb_mass( mrloop ) > 0 ) )
							{
								found[3] = true;
								mrloop_at_infall = mrloop;
							}
						}
					}
				}
				if ( mrloop_at_infall >= 3 )
				{
					if ( current_orbit.get_mass( mrloop_at_infall - 3 ) <= 0 )
					{
						#pragma omp atomic
						dropped++;
						if(debug_info)
						{
							std::cout << "Dropped orbit " << loop << endl;
						}
						if ( quiet == false )
							cout
									<< "WARNING: Dropped an orbit for insufficient pre-accretion lifetime."
									<< endl;
						continue;
					}
				}
				if ( current_orbit.get_mass( mrloop_at_infall )
						< mass_at_infall_cut )
				{
					#pragma omp atomic
					dropped++;
					if(debug_info)
					{
						std::cout << "Dropped orbit " << loop << endl;
					}
					if ( quiet == false )
						cout << "WARNING: Cut an orbit for low mass at infall."
								<< endl;
					continue;
				}
				if ( use_infalling == true )
				{
					//include halos with no robust infall time at a=1
					found[3] = true;
				}
				else
				{
					//include halos with no robust infall time as interlopers (preferred)
					if ( found[3] == false )
					{
						fmret = 1;
						found[3] = true;
					}
				}

				bool bad_result = false;

				if ( fmret < 0 )
				{
					brgastro::stripping_orbit temp_orbit_spline;
					if (  current_orbit.load_orbit_spline( mrloop_at_infall, &temp_orbit_spline, quiet )  )
					{
						fmret = 1;
						if ( (!debug_info) || (quiet==false) )
						{
							cout
									<< "WARNING: Could not load orbit " << loop << " for stripping calculation."
									<< endl;
						}
						bad_result = true;
					}
					else
					{
						bool record_this_orbit = record_all_orbits;
						if ( record_random_orbits )
						{
							if ( brgastro::drand<double>( 0., 1. ) < record_random_chance )
								record_this_orbit = true;
						}
						current_orbit.stripping_orbit_spline()->set_record_full_data(
								record_this_orbit );

						if ( debug_info )
						{
							cout
									<< "Trying to calculate stripping for orbit " << loop
									<< endl;
						}

						current_orbit.stripping_orbit_spline()->calc(!debug_info);
						if ( current_orbit.stripping_orbit_spline()->get_final_frac_m_ret(
								fmret ) )
						{
							bad_result = true;
							if(current_orbit.stripping_orbit_spline()->likely_disrupted())
								fmret = 0;
							else
								fmret = 1;
							if ( (debug_info) or (quiet==false) )
							{
								cout
										<< "WARNING: Could not calculate stripping for orbit " << loop
										<< endl;
							}
						}

						if ( debug_info )
							cout
									<< "Outputting stripping results for orbit " << loop
									<< endl;
						try {
							#pragma omp critical(recordorbits)
							if ( record_this_orbit )
							{
								std::vector< bool > output_params( 4, false );
								output_params[3] = true;
								current_orbit.stripping_orbit_spline()->set_satellite_output_parameters(
										output_params );
								output_params[0] = true;
								output_params[1] = true;
								output_params[3] = false;
								current_orbit.stripping_orbit_spline()->set_host_output_parameters(
										output_params );
								stringstream ss( "" );
								ss << output_file_name_base << output_file_count
										<< output_file_name_tail;
								string output_file_name = ss.str();
								ofstream out( output_file_name.c_str() );
								current_orbit.stripping_orbit_spline()->print_full_data(
										&out );
								out.close();
								out.clear();
								cout << " " << output_file_count << " " << loop << " "
										<< fmret << endl;
								output_file_count++;
							}
							#pragma omp critical(recordgeneral)
							if (record_general_orbit_data)
							{

								if ((current_orbit.stripping_orbit_spline()->likely_disrupted()) || (!current_orbit.stripping_orbit_spline()->bad_result())) {
									try
									{
										general_orbit_data_file << loop << "\t" << current_orbit.stripping_orbit_spline()->init_satellite_ptr()->mvir()*unitconv::kgtottMsun <<
												"\t" << current_orbit.stripping_orbit_spline()->final_host()->mvir()*unitconv::kgtottMsun << "\t" <<
												current_orbit.rf( t_now )*current_orbit.stripping_orbit_spline()->final_host()->rvir()*unitconv::mtokpc <<
												"\t" << current_orbit.rf( t_now ) << "\t" <<
												current_orbit.rfproj( t_now, 0 )*current_orbit.stripping_orbit_spline()->final_host()->rvir()*unitconv::mtokpc <<
												"\t" << current_orbit.rfproj( t_now, 0 ) << "\t" <<
												fabs(current_orbit.vf( t_now ))*current_orbit.stripping_orbit_spline()->final_host()->vvir()*unitconv::mpstokmps <<
												"\t" << fabs(current_orbit.vf( t_now )) << "\t" <<
												fabs(current_orbit.vfproj( t_now, 0 ))*current_orbit.stripping_orbit_spline()->final_host()->vvir()*unitconv::mpstokmps <<
												"\t" << fabs(current_orbit.vfproj( t_now, 0 )) << "\t" <<
												current_orbit.stripping_orbit_spline()->t_min_natural_value()*unitconv::stoGyr << "\t" <<
												current_orbit.stripping_orbit_spline()->last_infall_time()*unitconv::stoGyr << "\t" <<
												brgastro::zft(current_orbit.stripping_orbit_spline()->t_min_natural_value()) << "\t" <<
												brgastro::zft(current_orbit.stripping_orbit_spline()->last_infall_time()) << "\t" <<
												fmret << "\n";
										general_orbit_data_file.flush();
									}
									catch (const char * &e)
									{
										std::cerr << "WARNING: Could not output general orbit data for orbit #" << loop << ".\n"
												<< "Exception: " << e;
									}
									catch (const std::exception &e)
									{
										std::cerr << "WARNING: Could not output general orbit data for orbit #" << loop << ".\n"
												<< "Exception: " << e.what();
									}
								}
							}
						} catch (const std::exception &e) {
							std::cerr << "WARNING: Could not record orbit data for orbit " << loop << ".\n";
							bad_result = true;
						}

						if ( debug_info )
						{
							cout
									<< "Done outputting stripping results for orbit " << loop
									<< endl;
						}
					}
				}
				//add a count to the correct scale factor bin in the orbits R-V bin
				if ( (found[3] == true) && (!bad_result) )
				{
					if ( fmret >= 1 )
					{
						#pragma omp atomic
						interloper_count++;
						fmret = 1; // Just to ensure that if any mistakes slip through somehow, they're treated as interlopers
					}
					for ( int mrloop = 0; mrloop < nsteps; mrloop++ )
					{
						if ( ( fmret > mrbins[mrloop] )
								&& ( fmret <= mrbins[mrloop + 1] ) )
						{
							mr = mrloop;
							found[4] = true;
							break;
						}
					}
					if ( found[4] == true )
					{
						#pragma omp atomic
						pdfs[v][r][mr]++;
					}
					else
					{
						#pragma omp atomic
						dropped++;
						if ( quiet == false )
							cout << "WARNING: No bin found for orbit." << endl;
						continue;
					}
				}
				else
				{
					#pragma omp atomic
					dropped++;
					if ( quiet == false )
						cout << "WARNING: No infall time computed for orbit."
								<< endl;
					continue;
				}
			}
			else
			{
				#pragma omp atomic
				dropped++;
				int t_now = current_orbit.get_nsteps() - 1;
				if ( quiet == false )
					cout << "WARNING: No R-V bin(s) found for orbit. R="
							<< current_orbit.rfproj( t_now, 0 ) << " V="
							<< abs( current_orbit.vfproj( t_now, 0 ) ) << endl;
				continue;
			}

			try {
				current_orbit.unload_orbit_spline();
			} catch (const std::exception &e) {
				std::cerr << "WARNING: Could not unload orbit spline!\n";
			}

		}

		cout << "Dropped " << dropped << "/" << Orbits.size() << " orbits."
				<< endl;
	} while((Orbits.size() >= orbit_chunk_size) && (orbit_chunk_size>0));

	if(record_general_orbit_data) general_orbit_data_file.close();

	cout << "Writing pdf table...";
	cout.flush();

	for ( int vloop = 0; vloop < Nbinsv; vloop++ )
		for ( int rloop = 0; rloop < Nbinsr; rloop++ )
			for ( int aloop = 0; aloop < nsteps; aloop++ )
				outfile << " " << pdfs[vloop][rloop][aloop];

	outfile.flush();
	outfile.close();

	cout << "Done." << endl;

	//free pdfs memory:
	for ( int oloop = 0; oloop < Nbinsv; oloop++ )
	{
		for ( int loop = 0; loop < Nbinsr; loop++ )
			delete[] pdfs[oloop][loop];
	}
	for ( int loop = 0; loop < Nbinsv; loop++ )
	{
		delete[] pdfs[loop];
	}

	delete[] pdfs;

	return 0;
}

