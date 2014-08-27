/**********************************************************************\
  @file SALTSA_file_functions.h

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

// body file: SALTSA_file_functions.cpp

#ifndef SALTSA_FILE_FUNCTIONS_H_
#define SALTSA_FILE_FUNCTIONS_H_

#include <iostream>
#include <string>
#include <vector>

#include "SALTSA_global.h"

namespace SALTSA {

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

} // namespace SALTSA

#endif /* SALTSA_FILE_FUNCTIONS_H_ */
