#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <new>
#include <fstream>
#include <string>
#include "SALTSA_misc_functions.h"
#include "SALTSA_calculus.hpp"

using namespace std;

/** Class method definitions **/
#if (1)
// SALTSA::phase function implementations
#if (1)
/**
 *
 * @param init_x
 * @param init_y
 * @param init_z
 * @param init_vx
 * @param init_vy
 * @param init_vz
 * @param init_t
 */
SALTSA::phase::phase( double init_x, double init_y,
double init_z,
double init_vx, double init_vy, double init_vz,
double init_t )
{
	x = init_x;
	y = init_y;
	z = init_z;
	vx = init_vx;
	vy = init_vy;
	vz = init_vz;
	t = init_t;
}

/**
 *
 * @param init_x
 * @param init_y
 * @param init_z
 * @param init_vx
 * @param init_vy
 * @param init_vz
 * @param init_t
 * @return
 */
const int SALTSA::phase::set_phase( double init_x, double init_y,
double init_z,
double init_vx, double init_vy, double init_vz,
double init_t )
{
	x = init_x;
	y = init_y;
	z = init_z;
	vx = init_vx;
	vy = init_vy;
	vz = init_vz;
	t = init_t;

	return 0;
}

#endif // end SALTSA::phase function implementations

#endif // end class function definitions

/** Global function implementations **/
#if (1)
/**
 *
 * @param value
 * @param epsilon
 * @return
 */
const int SALTSA::round_int( const double value, const double epsilon )
{

	if ( value < 0.0 )
		return -SALTSA::round_int( -value, epsilon );

	double ipart;
	std::modf( value, &ipart );

	// If 'value' is exctly halfway between two integers
	if ( fabs( value - ( ipart + 0.5 ) ) < epsilon )
	{
		// If 'ipart' is even then return 'ipart'
		if ( std::fmod( ipart, 2.0 ) < epsilon )
			return (int)ipart;

		// Else return the nearest even integer
		return (int)ceil( ipart + 0.5 );
	}

	// Otherwise use the usual round to closest
	// (Either symmetric half-up or half-down will do0
	return (int)floor( value + 0.5 );
}

/**
 *
 * @param out_stream
 * @param num_columns
 * @param num_rows
 * @param header
 * @param data
 * @param skip_header
 * @param silent
 * @return
 */
const int SALTSA::print_table( std::ostream & out_stream,
		const int num_columns, const int num_rows,
		const std::vector< std::string > & header,
		const std::vector< std::vector< std::string > > &data,
		const bool skip_header, const bool silent )
{
	vector< int > width;
	if ( int errcode = SALTSA::make_array( width, num_columns ) )
		return errcode + LOWER_LEVEL_ERROR;

	try
	{

		// First, we loop through to get the maximum width of each column
		// Check the header first
		for ( int c = 0; c < num_columns; c++ )
		{
			if ( (int)header[c].length() > width[c] )
			{
				width[c] = header[c].length();
			}
		} // for( int c = 0; c < num_columns; c++ )

		// Now loop through the data
		for ( int i = 0; i < num_rows; i++ )
		{
			for ( int c = 0; c < num_columns; c++ )
			{
				if ( (int)data[c][i].length() > width[c] )
				{
					width[c] = data[c][i].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		} // for( int i = 0; i < num_rows; i++ ) (testing width)

		// Increase all widths by 1 to ensure spacing
		for ( int c = 0; c < num_columns; c++ )
			width[c] += 1;

		// Output the header
		if ( !skip_header )
			for ( int c = 0; c < num_columns; c++ )
				out_stream << setfill( ' ' ) << setw( width[c] ) << header[c];

		out_stream << std::endl;

		// Output the data
		for ( int i = 0; i < num_rows; i++ )
		{
			for ( int c = 0; c < num_columns; c++ )
			{
				out_stream << setfill( ' ' ) << setw( width[c] ) << data[c][i];
			}
			out_stream << std::endl;
		}

	}
	catch ( ... )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Could not print table. Check that the data is properly formatted\n"
					<< "at least num_columns length for header and first index of data, and at\n"
					<< "least num_rows length for all vectors contained within data.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	return 0;
}

/**
 *
 * @param out_stream
 * @param num_columns
 * @param num_rows
 * @param data
 * @param silent
 * @return
 */
const int SALTSA::print_table( std::ostream & out_stream,
		const int num_columns, const int num_rows,
		const std::vector< std::vector< std::string > > & data,
		const bool silent )
{
	vector< int > width;
	if ( int errcode = SALTSA::make_array( width, num_columns ) )
		return errcode + LOWER_LEVEL_ERROR;

	try
	{

		// First, we loop through to get the maximum width of each column
		// Loop through the data
		for ( int i = 0; i < num_rows; i++ )
		{
			for ( int c = 0; c < num_columns; c++ )
			{
				if ( (int)data[c][i].length() > width[c] )
				{
					width[c] = data[c][i].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		} // for( int i = 0; i < num_rows; i++ ) (testing width)

		// Increase all widths by 1 to ensure spacing
		for ( int c = 0; c < num_columns; c++ )
			width[c] += 1;

		out_stream << std::endl;

		// Output the data
		for ( int i = 0; i < num_rows; i++ )
		{
			for ( int c = 0; c < num_columns; c++ )
			{
				out_stream << setfill( ' ' ) << setw( width[c] ) << data[c][i];
			}
			out_stream << std::endl;
		}

	}
	catch ( ... )
	{
		if ( !silent )
			std::cerr
					<< "ERROR: Could not print table. Check that the data is properly formatted\n"
					<< "at least num_columns length for first index of data, and at\n"
					<< "least num_rows length for all vectors contained within data.\n";
		return INVALID_ARGUMENTS_ERROR;
	}

	return 0;
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
/**
 *
 * @param table_file_name
 * @param table_data
 * @param silent
 * @return
 */
const int SALTSA::load_table( const std::string & table_file_name, std::vector<std::vector<double> > & table_data,
		const bool silent)
{
	std::ifstream fi;

	// Open the file
	if(SALTSA::open_file(fi, table_file_name))
		throw std::runtime_error("Cannot open file in load_table.");

	// Trim the header
	if(trim_comments_all_at_top(fi))
		throw std::runtime_error("Invalid table format in load_table");

	// Clear the output vector
	table_data.resize(0);

	std::string line_data;
	std::istringstream line_data_stream;
	while ( !fi.eof() )
	{
		std::vector<double> temp_vector(0);

		getline(fi, line_data);
		line_data_stream.clear();
		line_data_stream.str(line_data);

	    // Split the line on whitespace
		do
	    {
	        double value;
	        line_data_stream >> value;
	        temp_vector.push_back(value);
	    } while (line_data_stream);

	    table_data.push_back(temp_vector);
	}

	fi.close();
	fi.clear();

	return 0;
}
/**
 *
 * @param table_file_name
 * @param table_data
 * @param header
 * @param silent
 * @return
 */
const int SALTSA::load_table( const std::string & table_file_name,
		std::vector<std::vector<double> > & table_data,
		std::vector<std::string> & header, const bool silent)
{
	std::ifstream fi;
	std::string line_data;
	std::istringstream line_data_stream;

	// Open the file
	if(SALTSA::open_file(fi, table_file_name))
		throw std::runtime_error("Cannot open file in load_table.");

	// Clear the output vectors
	header.resize(0);
	table_data.resize(0);

	// Get the header line and store it
	getline(fi, line_data);
	line_data_stream.clear();
	line_data_stream.str(line_data);

    // Split the line on whitespace
	do
    {
        std::string value;
        line_data_stream >> value;
        header.push_back(value);
    } while (line_data_stream);

	// Check if the first entry in the header is a comment marker. If so, remove it
	if(header.at(0) == "#")
	{
		header.erase(header.begin());
	}
	else
	{
		// Check if the first character of the first entry is a comment marker. Again,
		// delete it if so (but not the whole entry!)
		if(header.at(0).at(0) == (char)'#')
		{
			header.at(0).erase(header.at(0).begin());
		}
	}

	// Trim any remaining comments
	if(trim_comments_all_at_top(fi))
		throw std::runtime_error("Invalid table format in load_table");

	// Load in the data now
	while ( !fi.eof() )
	{
		std::vector<double> temp_vector(0);

		getline(fi, line_data);
		if(line_data.size()==0) break;
		line_data_stream.clear();
		line_data_stream.str(line_data);

	    // Split the line on whitespace
		do
	    {
	        double value;
	        line_data_stream >> value;
	        temp_vector.push_back(value);
	    } while (line_data_stream);

	    table_data.push_back(temp_vector);
	}

	fi.close();
	fi.clear();

	return 0;
}
/**
 *
 * @param table_file_name
 * @param header_links
 * @param case_sensitive
 * @param silent
 * @return
 */
const int SALTSA::load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<double>* > > & header_links,
		const bool case_sensitive, const bool silent)
{
	// First, load in the table
	std::vector<std::vector<double> > table_data;
	std::vector<std::string> header;

	SALTSA::load_table( table_file_name, table_data, header, silent);

	// Now, loop through each key and search for it in the header.
	for(unsigned int i = 0; i < header_links.size(); i++)
	{
		std::string test_string_key( header_links.at(i).first );

		if(!case_sensitive)
		{
			std::transform(test_string_key.begin(), test_string_key.end(),
					test_string_key.begin(), ::tolower);
		}

		bool key_found = false;

		for(unsigned int j = 0; j < header.size(); j++)
		{
			std::string test_string_header(header.at(j));

			if(!case_sensitive)
			{
				std::transform(test_string_header.begin(), test_string_header.end(),
						test_string_header.begin(), ::tolower);
			}

			if(test_string_key==test_string_header)
			{
				// Found it, now load in the data
				key_found = true;

				std::vector<double> column_data(table_data.size());
				try
				{
					for(unsigned int k = 0; k < table_data.size(); k++)
					{
						column_data.at(k) = table_data.at(k).at(j);
					}
					*(header_links.at(i).second) = column_data;
				}
				catch (std::exception &e)
				{
					throw std::runtime_error("Improperly formatted data table in load_table_columns.");
				}

				break;
			}
		}

		if(!key_found)
		{
			throw std::runtime_error("Key not found in header in load_table_columns.");
		}
	}

	return 0;
}

/**
 *
 * @param stream
 * @param name
 * @param silent
 * @return
 */
const int SALTSA::open_file( std::ofstream & stream, const std::string name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		if ( !silent )
			cerr << "ERROR: Could not open file " << name << "!\n";
		return FILE_ACCESS_ERROR;
	}
	return 0;
}
/**
 *
 * @param stream
 * @param name
 * @param silent
 * @return
 */
const int SALTSA::open_file( std::ifstream & stream, const std::string name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		if ( !silent )
			cerr << "ERROR: Could not open file " << name << "!\n";
		return FILE_ACCESS_ERROR;
	}
	return 0;
}
/**
 *
 * @param stream
 * @param name
 * @param silent
 * @return
 */
const int SALTSA::open_file( std::fstream & stream, const std::string name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		if ( !silent )
			cerr << "ERROR: Could not open file " << name << "!\n";
		return FILE_ACCESS_ERROR;
	}
	return 0;
}

/**
 *
 * @param stream
 * @param silent
 * @return
 */
const int SALTSA::trim_comments_one_line( std::ifstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
	else
	{
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

/**
 *
 * @param stream
 * @param silent
 * @return
 */
const int SALTSA::trim_comments_one_line( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
	else
	{
		return UNSPECIFIED_ERROR;
	}
	return 0;
}

/**
 *
 * @param stream
 * @param silent
 * @return
 */
const int SALTSA::trim_comments_all_at_top( std::ifstream & stream,
		const bool silent )
{
	std::string file_data;
	while ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
		{
			getline( stream, file_data );
		}
		else
		{
			return 0;
		}
	}
	return UNSPECIFIED_ERROR; // We reached the end before we ran out of comments
}

/**
 *
 * @param stream
 * @param silent
 * @return
 */
const int SALTSA::trim_comments_all_at_top( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	while ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
		{
			getline( stream, file_data );
		}
		else
		{
			return 0;
		}
	}
	return UNSPECIFIED_ERROR; // We reached the end before we ran out of comments
}

/**
 *
 * @param mean
 * @param stddev
 * @return
 */
const double SALTSA::Gaus_rand( double mean, double stddev )
{
	double x1, x2, w;

	if ( stddev <= 0 )
		return mean;

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( mean + x1 * w * stddev );

} // double Gaus_rand(double mean, double stddev)

/**
 *
 * @param mean
 * @param stddev
 * @return
 */
const double SALTSA::log10Gaus_rand( double mean, double stddev )
{
	double x1, x2, w, fact;

	if ( stddev <= 0 )
		return mean;

	stddev *= log( 10 ); // Converts dex to natural log
	fact = exp( -std::pow( stddev, 2 ) / 2 );

	do
	{
		x1 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		x2 = 2.0 * ( ( (double)rand() ) / RAND_MAX ) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( ( -2.0 * log( w ) ) / w );
	return ( mean * fact * exp( x1 * w * stddev ) );
} // double lnGaus_rand(double mean, double stddev)

/**
 *
 * @param lambda
 * @return
 */
const int SALTSA::Pois_rand( double lambda )
{
	double L, p;
	int k;
	L = exp( -lambda );
	k = 0;
	p = 1;
	do
	{
		k++;
		p *= drand( 0, 1 );
	} while ( p > L );

	return ( k - 1 );
} // int Pois_rand(double lambda)

#endif // end global function implementations
