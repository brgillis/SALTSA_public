/**********************************************************************\
  @file file_functions.h

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

// body file: brg/file_functions.cpp

#ifndef _BRG_FILE_FUNCTIONS_H_INCLUDED_
#define _BRG_FILE_FUNCTIONS_H_INCLUDED_

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "global.h"

#include "brg/utility.hpp"

namespace brgastro
{

// Typedefs
#if(1)

struct header {
    typedef std::vector< std::string > type;
};

template <typename T>
struct table {
    typedef std::vector< std::vector< T > > type;
};

template <typename T>
struct table_map {
    typedef std::map< std::string, std::vector< T > > type;
};

#endif

/** Global function declarations **/
#if (1)

// Functions to open a file and check that it's been opened successfully. An exception will be thrown
// if the file isn't opened successfully.
void open_file( std::ofstream & stream, const std::string & name,
		const bool silent = false );
void open_file( std::ifstream & stream, const std::string & name,
		const bool silent = false );
void open_file( std::fstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::ofstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::ifstream & stream, const std::string & name,
		const bool silent = false );
void open_bin_file( std::fstream & stream, const std::string & name,
		const bool silent = false );

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
void print_table( std::ostream & out_stream,
		const table<std::string>::type & data,
		const header::type & header = header::type(),
		const bool silent = false );

// Some templates to coerce any type of data to be printed out
template<typename T>
void print_table( std::ostream & out_stream,
		const typename table<T>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::stringstream ss;
	table<std::string>::type string_data(0);

	make_array2d(string_data, data);
	for( size_t i=0; i<data.size(); ++i )
	{
		for( size_t j=0; i<data[i].size(); ++j )
		{
			ss.str("");
			ss << data[i][j];
			string_data[i].at(j) = ss.str();
		}
	}
	print_table(out_stream, string_data, header, silent);
}

// And to allow us to print to a file name instead of a stream
inline void print_table( const std::string & file_name,
		const table<std::string>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::ofstream fo;
	open_file(fo,file_name,silent);

	print_table(fo,data,header,silent);
}
template<typename T>
void print_table( const std::string & file_name,
		const typename table<T>::type & data,
		const header::type & header = header::type(),
		const bool silent = false )
{
	std::ofstream of;
	open_file(of,file_name,silent);

	print_table<T>(of,data,header,silent);
}

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
table<std::string>::type load_table( std::istream & table_file_name,
		const bool silent=false);

template<typename T>
typename table<T>::type load_table( std::istream & table_file_name,
		const bool silent=false)
{
	std::stringstream ss;

	table<std::string>::type string_data(load_table( table_file_name, silent ));
	typename table<T>::type data;

	make_array2d(data,string_data);

	for(size_t i=0; i<string_data.size(); ++i)
	{
		for(size_t j=0; j<string_data[i].size(); ++j)
		{
			ss.clear();
			ss.str(string_data[i][j]);
			ss >> data[i].at(j);
		}
	}
	return data;
}

// And to allow us to load from a file name instead of a stream
inline table<std::string>::type load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table(fi,silent);
}
template<typename T>
inline typename table<T>::type load_table( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table<T>(fi,silent);
}

// Load a table's header as a vector of strings
header::type load_header( std::istream & table_stream,
		const bool silent=false);

// And to allow us to load from a file name instead of a stream
inline header::type load_header( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_header(fi,silent);
}

// Merge a header and data table into a map
template<typename T>
typename table_map<T>::type make_table_map(
		const typename table<T>::type & data,
		const header::type & header,
		const bool silent=false)
{
	typename table_map<T>::type result;

	size_t h_size = header.size();
	size_t d_size = data.size();

	const header::type *header_to_use = &header;
	header::type new_header;
	const typename table<T>::type *data_to_use = &data;
	typename table<T>::type new_data;

	// First, check if the header and data aren't the same size
	if(h_size<d_size)
	{
		// We'll pad the header with default values
		new_header = header;
		new_header.resize(d_size);
		for(size_t i=h_size; i<d_size; ++i)
		{
			std::stringstream ss("col_");
			ss << i+1;
			new_header[i] = ss.str();
		}
		header_to_use = &new_header;
		h_size = d_size;
	}
	else if(d_size<h_size)
	{
		// We'll pad the data with default values
		new_data = data;
		new_data.resize(h_size);
		if(d_size>0)
		{
			// If we have some data, match the size of it in new columns
			for(size_t i=d_size; i<h_size; ++i)
			{
				make_array1d(new_data[i],new_data[0].size());
			}
		}
		data_to_use = &new_data;
		d_size = h_size;
	}

	for(size_t i=0; i<h_size; ++i)
	{
		result[(*header_to_use)[i]] = (*data_to_use)[i];
	}
	return result;
}
inline table_map<std::string>::type make_table_map(
		const table<std::string>::type & data,
		const header::type & header,
		const bool silent=false)
{
	return make_table_map<std::string>(data,header,silent);
}

// Directly load a map of the data
inline table_map<std::string>::type load_table_map( std::istream & fi,
		const bool silent=false)
{
	header::type header = load_header(fi,silent);
	table<std::string>::type data = load_table(fi,silent);
	return make_table_map(data,header,silent);
}

template<typename T>
typename table_map<T>::type load_table_map( std::istream & fi,
		const bool silent=false)
{
	header::type header = load_header(fi,silent);
	typename table<T>::type data = load_table<T>(fi,silent);
	return make_table_map<T>(data,header,silent);
}

// And to allow us to load from a file name instead of a stream
inline table_map<std::string>::type load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_map(fi,silent);
}
template<typename T>
typename table_map<T>::type load_table_map( const std::string & file_name,
		const bool silent = false )
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	return load_table_map<T>(fi,silent);
}

inline void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<std::string>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	table_map<std::string>::type table_map = load_table_map(fi);

	for(std::map< std::string, std::vector<std::string>* >::iterator it=column_map.begin();
			it!=column_map.end(); ++it)
	{
		*(it->second) = table_map[it->first];

		// Check that we found it
		if(!silent)
		{
			if(it->second->size()==0)
			{
				std::cerr << "WARNING: Column " << it->first << " not found in table.\n";
			}
		}

	}
}

template<typename T>
void load_table_columns( std::istream & fi,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	typename table_map<T>::type table_map = load_table_map<T>(fi);

	for(typename std::map< std::string, std::vector<T>* >::iterator it=column_map.begin();
			it!=column_map.end(); ++it)
	{
		*(it->second) = table_map[it->first];

		// Check that we found it
		if(!silent)
		{
			if(it->second->size()==0)
			{
				std::cerr << "WARNING: Column " << it->first << " not found in table.\n";
			}
		}

	}
}

// And to allow us to load from a file name instead of a stream
inline void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<std::string>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	load_table_columns(fi,column_map,case_sensitive,silent);
}
template<typename T>
void load_table_columns( const std::string & file_name,
		std::map< std::string, std::vector<T>* > & column_map,
		const bool case_sensitive=false, const bool silent=false)
{
	std::ifstream fi;
	open_file(fi,file_name,silent);

	load_table_columns<T>(fi,column_map,case_sensitive,silent);
}

// Functions to get rid of comments lines (those starting with #) from open fstreams and ifstreams. The "one_line" versions will
// trim the next line if it's commented out, and the "all_at_top" versions will trim all commented lines from the current
// position until they find a line which isn't a comment.
// The "one_line" functions will return 1 if the file is already at the end, and 0 otherwise.
// The "all_at_top" functions will return 1 if they reach the end of the file before the run out of comments, and 0 otherwise.
void trim_comments_one_line( std::istream & stream, const bool silent =
		false );
void trim_comments_one_line( std::fstream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::istream & stream, const bool silent =
		false );
void trim_comments_all_at_top( std::fstream & stream, const bool silent =
		false );

// Utility functions

// Splits a string into a vector of "word" strings on whitespace
std::vector< std::string > split_on_whitespace( const std::string & sentence );
header::type convert_to_header( const std::string & line );

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
