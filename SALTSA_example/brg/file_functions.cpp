/**********************************************************************\
  @file file_functions.cpp

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

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>

#include "brg/global.h"

#include "brg/utility.hpp"

#include "brg/file_functions.h"

using namespace std;

/** Global function implementations **/
#if (1)
void brgastro::_print_table( std::ostream & out_stream,
		const std::vector< std::vector< std::string > > &data,
		const std::vector< std::string > & header,
		const bool silent )
{
	unsigned int num_columns = data.size();
	unsigned int num_rows = data.at(0).size();
	vector< int > width(num_columns,0);

	const bool skip_header = (header.size()==0);

	try
	{
		// First, we loop through to get the maximum width of each column
		// Check the header first
		if(!skip_header)
		{
			for ( unsigned int c = 0; c < num_columns; c++ )
			{
				if ( (int)header[c].length() > width[c] )
				{
					width[c] = header[c].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		}

		// Now loop through the data
		for ( unsigned int i = 0; i < num_rows; i++ )
		{
			for ( unsigned int c = 0; c < num_columns; c++ )
			{
				if ( (int)data[c].at(i).length() > width[c] )
				{
					width[c] = data[c][i].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		} // for( int i = 0; i < num_rows; i++ ) (testing width)

		// Increase all widths by 1 to ensure spacing
		for ( unsigned int c = 0; c < num_columns; c++ )
			width[c] += 1;

		// Output the header
		if ( !skip_header )
			for ( unsigned int c = 0; c < num_columns; c++ )
				out_stream << setfill( ' ' ) << setw( width[c] ) << header[c];

		out_stream << std::endl;

		// Output the data
		for ( unsigned int i = 0; i < num_rows; i++ )
		{
			for ( unsigned int c = 0; c < num_columns; c++ )
			{
				out_stream << setfill( ' ' ) << setw( width[c] ) << data[c][i];
			}
			out_stream << std::endl;
		}

	}
	catch ( const std::out_of_range &e )
	{
		throw std::runtime_error((std::string)"ERROR: Could not print table. Check that the data is properly formatted\n"
					+ "at least num_columns length for header and first index of data, and at\n"
					+ "least num_rows length for all vectors contained within data.\n");
	}
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
std::vector<std::vector<std::string> > brgastro::_load_table( const std::string & table_file_name, const bool silent)
{
	std::vector<std::vector<std::string> > table_data;
	std::ifstream fi;

	// Open the file
	open_file(fi, table_file_name);

	// Trim the header
	trim_comments_all_at_top(fi);

	// Clear the output vector
	table_data.resize(0);

	std::string line_data;
	std::istringstream line_data_stream;
	while ( !fi.eof() )
	{
		std::vector<std::string> temp_vector(0);

		getline(fi, line_data);
		line_data_stream.clear();
		line_data_stream.str(line_data);

	    // Split the line on whitespace
		do
	    {
	        std::string value;
	        line_data_stream >> value;
	        temp_vector.push_back(value);
	    } while (line_data_stream);

	    table_data.push_back(temp_vector);
	}

	return table_data;
}
void brgastro::_load_table_and_header( const std::string & table_file_name,
		std::vector<std::vector<std::string> > & table_data,
		std::vector<std::string> & header, const bool silent)
{
	std::ifstream fi;
	std::string line_data;
	std::istringstream line_data_stream;

	// Open the file
	brgastro::open_file(fi, table_file_name);

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
	trim_comments_all_at_top(fi);

	// Load in the data now
	while ( !fi.eof() )
	{
		std::vector<std::string> temp_vector(0);

		getline(fi, line_data);
		if(line_data.size()==0) break;
		line_data_stream.clear();
		line_data_stream.str(line_data);

	    // Split the line on whitespace
		do
	    {
	        std::string value;
	        line_data_stream >> value;
	        temp_vector.push_back(value);
	    } while (line_data_stream);

	    table_data.push_back(temp_vector);
	}
}
void brgastro::_load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<std::string>* > > & header_links,
		const bool case_sensitive, const bool silent)
{
	// First, load in the table
	std::vector<std::vector<std::string> > table_data;
	std::vector<std::string> header;

	brgastro::_load_table_and_header( table_file_name, table_data, header, silent);

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

				std::vector<std::string> column_data(table_data.size());
				try
				{
					for(unsigned int k = 0; k < table_data.size(); k++)
					{
						column_data.at(k) = table_data.at(k).at(j);
					}
					*(header_links.at(i).second) = column_data;
				}
				catch (const std::out_of_range &e)
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

	return;
}

void brgastro::open_file( std::ofstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_file( std::ifstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_file( std::fstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str() );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_bin_file( std::ofstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), ios::out | ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_bin_file( std::ifstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), ios::in | ios::binary );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}
void brgastro::open_bin_file( std::fstream & stream, const std::string & name,
		const bool silent )
{
	stream.close();
	stream.clear();
	stream.open( name.c_str(), ios::out | ios::in | ios::binary  );
	if ( !stream )
	{
		throw std::runtime_error("ERROR: Could not open file " + name + "!\n");
	}
}

void brgastro::trim_comments_one_line( std::ifstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_one_line( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( !stream.eof() )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_all_at_top( std::ifstream & stream,
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
			return;
		}
	}
}

void brgastro::trim_comments_all_at_top( std::fstream & stream,
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
			return;
		}
	}
}

#endif // end global function implementations
