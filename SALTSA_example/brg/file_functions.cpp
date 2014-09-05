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
#include "brg/vector/manipulations.hpp"

#include "file_functions.h"

using namespace std;

/** Global function implementations **/
#if (1)
void brgastro::print_table( std::ostream & out_stream,
		const std::vector< std::vector< std::string > > &data,
		const std::vector< std::string > & header,
		const bool silent )
{
	size_t num_columns = data.size();
	size_t num_rows = data.at(0).size();
	vector< int > width(num_columns,0);

	const bool skip_header = (header.size()==0);

	try
	{
		// First, we loop through to get the maximum width of each column
		// Check the header first
		if(!skip_header)
		{
			for ( size_t c = 0; c < num_columns; c++ )
			{
				if ( (int)header[c].length() > width[c] )
				{
					width[c] = header[c].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		}

		// Now loop through the data
		for ( size_t i = 0; i < num_rows; i++ )
		{
			for ( size_t c = 0; c < num_columns; c++ )
			{
				if ( (int)data[c].at(i).length() > width[c] )
				{
					width[c] = data[c][i].length();
				}
			} // for( int c = 0; c < num_columns; c++ )
		} // for( int i = 0; i < num_rows; i++ ) (testing width)

		// Increase all widths by 1 to ensure spacing
		for ( size_t c = 0; c < num_columns; c++ )
			width[c] += 1;

		// Output the header
		if ( !skip_header )
			for ( size_t c = 0; c < num_columns; c++ )
				out_stream << setfill( ' ' ) << setw( width[c] ) << header[c];

		out_stream << std::endl;

		// Output the data
		for ( size_t i = 0; i < num_rows; i++ )
		{
			for ( size_t c = 0; c < num_columns; c++ )
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
std::vector<std::vector<std::string> > brgastro::load_table( std::istream & fi,
		const bool silent)
{
	std::vector<std::vector<std::string> > table_data;

	// Trim the header
	trim_comments_all_at_top(fi);

	// Clear the output vector
	table_data.resize(0);

	std::string line_data;
	std::istringstream line_data_stream;
	while ( getline(fi, line_data) )
	{
		std::vector<std::string> temp_vector(0);

		line_data_stream.clear();
		line_data_stream.str(line_data);

	    // Split the line on whitespace
        std::string value;
		while (line_data_stream >> value)
	    {
	        temp_vector.push_back(value);
	    }

	    table_data.push_back(temp_vector);
	}

	return transpose(table_data);
}

brgastro::header::type brgastro::load_header( std::istream & table_stream,
		const bool silent)
{
	std::string temp_line;
	std::vector<header::type> possible_headers;

	// Get all comment lines at the top of the file
	while ( table_stream )
	{
		if ( table_stream.peek() == (int)( *"#" ) )
		{
			getline( table_stream, temp_line );

			header::type temp_header = convert_to_header(temp_line);

			if(temp_header.size() > 0)
				possible_headers.push_back(temp_header);
		}
		else
		{
			break;
		}
	}

	if(possible_headers.size()==1) return possible_headers[0];
	if(possible_headers.size()==0) return header::type();

	// If we get here, more than one line is a possible header candidate. Our next step is to
	// go to the data and count the columns in the first line. If only one possible header has
	// the right length, we know that's the one.

	unsigned short int n_cols = 0;
	do
	{
		getline( table_stream, temp_line );
		std::string junk;
		std::istringstream line_data_stream(temp_line);
		while (line_data_stream >> junk)
		{
			++n_cols;
		}
	} while(n_cols==0 && table_stream);

	if(n_cols==0) // If we can't find any data
	{
		std::cerr << "ERROR: Header line ambiguous; returning null header.\n";
		return header::type();
	}

	// Search through the possible headers, and see if we find exactly one with the right size
	unsigned short int num_right_size = 0;
	size_t i_best = 0;
	for(size_t i=0; i<possible_headers.size(); ++i)
	{
		if(possible_headers[i].size()==n_cols)
		{
			++num_right_size;
			i_best = i;
		}
	}

	if(num_right_size != 1) // If multiple or zero lines are the right size
	{
		std::cerr << "ERROR: Header line ambiguous; returning null header.\n";
		return header::type();
	}

	return possible_headers[i_best];
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

void brgastro::trim_comments_one_line( std::istream & stream,
		const bool silent )
{
	std::string file_data;
	if ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_one_line( std::fstream & stream,
		const bool silent )
{
	std::string file_data;
	if ( stream )
	{
		if ( stream.peek() == (int)( *"#" ) )
			getline( stream, file_data );
	}
}

void brgastro::trim_comments_all_at_top( std::istream & stream,
		const bool silent )
{
	std::string file_data;
	while ( stream )
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
	while ( stream )
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

std::vector< std::string > brgastro::split_on_whitespace( const std::string & sentence )
{
	std::vector< std::string > result;
	std::istringstream sentence_data_stream(sentence);

	std::string word;
	while (sentence_data_stream >> word)
	{
		result.push_back(word);
	}

	return result;
}
brgastro::header::type brgastro::convert_to_header( const std::string & line )
{
	header::type result;
	std::istringstream line_data_stream(line);

	std::string word;

	// Get rid of first word if it's the comment indicator
	if ( line_data_stream.peek() == (int)( *"#" ) )
	{
		line_data_stream >> word;
	}

	while (line_data_stream >> word)
	{

		result.push_back(word);
	}

	return result;
}

#endif // end global function implementations
