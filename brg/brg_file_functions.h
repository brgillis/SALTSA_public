/**********************************************************************\
brg_file_functions.h
 -----------

 If this header is used, the source file brg_file_functions.cpp must be
 included and compiled with the project.

 This file includes various functions having to do with file reading and
 writing.

 Everything in this file is declared in the namespace brgastro.

 \**********************************************************************/

#ifndef __BRG_FILE_FUNCTIONS_H_INCLUDED__
#define __BRG_FILE_FUNCTIONS_H_INCLUDED__

#include "brg_global.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>

#include "brg_units.h"
#include "brg_misc_functions.hpp"

namespace brgastro
{

/** Global function declarations **/
#if (1)

// Prints a formatted table in the passed stream. header is a vector of strings representing the labels for each column,
// and data is a 2-d vector of the data to be printed, in the format data[c][r], where c is the column index and r is the row index.
// Prints an error message and returns 1 if an error arises (for instance, vectors are too short).
const int print_table( std::ostream & out_stream, const int num_columns,
		const int num_rows, const std::vector< std::string > & header,
		const std::vector< std::vector< std::string > > & data,
		const bool skip_header = false, const bool silent = false );
const int print_table( std::ostream & out_stream, const int num_columns,
		const int num_rows,
		const std::vector< std::vector< std::string > > & data,
		const bool silent = false );

// Load table, either loading in the entire table, or only loading in certain columns into pointed-to
// variables, found by matching header entries to the strings passed
const int load_table( const std::string & table_file_name, std::vector<std::vector<double> > & table_data,
		const bool silent=false);
const int load_table( const std::string & table_file_name, std::vector<std::vector<double> > & table_data,
		std::vector<std::string> & header, const bool silent=false);
const int load_table_columns( const std::string & table_file_name,
		std::vector< std::pair< std::string, std::vector<double>* > > & header_links,
		const bool case_sensitive=false, const bool silent=false);

// Functions to open a file and check that it's been opened successfully. An error message will be printed if it can't be opened,
// and the functions will return 1. Otherwise they'll return 0.
const int open_file( std::ofstream & stream, const std::string name,
		const bool silent = false );
const int open_file( std::ifstream & stream, const std::string name,
		const bool silent = false );
const int open_file( std::fstream & stream, const std::string name,
		const bool silent = false );
const int open_bin_file( std::ofstream & stream, const std::string name,
		const bool silent = false );
const int open_bin_file( std::ifstream & stream, const std::string name,
		const bool silent = false );
const int open_bin_file( std::fstream & stream, const std::string name,
		const bool silent = false );

// Functions to get rid of comments lines (those starting with #) from open fstreams and ifstreams. The "one_line" versions will
// trim the next line if it's commented out, and the "all_at_top" versions will trim all commented lines from the current
// position until they find a line which isn't a comment.
// The "one_line" functions will return 1 if the file is already at the end, and 0 otherwise.
// The "all_at_top" functions will return 1 if they reach the end of the file before the run out of comments, and 0 otherwise.
const int trim_comments_one_line( std::ifstream & stream, const bool silent =
		false );
const int trim_comments_one_line( std::fstream & stream, const bool silent =
		false );
const int trim_comments_all_at_top( std::ifstream & stream, const bool silent =
		false );
const int trim_comments_all_at_top( std::fstream & stream, const bool silent =
		false );

#endif // End global function declarations

}

#endif // __BRG_FILE_FUNCTIONS_H_INCLUDED__
