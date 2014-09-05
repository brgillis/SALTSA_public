/**********************************************************************\
 @file manipulations.hpp
 ------------------

 Functions to manipulate vectors in various ways.

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

#ifndef _BRG_MANIPULATIONS_HPP_INCLUDED_
#define _BRG_MANIPULATIONS_HPP_INCLUDED_

#include <vector>

#include "brg/global.h"

#include "brg/utility.hpp"

namespace brgastro {

template<typename T>
std::vector< std::vector<T> > transpose(const std::vector< std::vector<T> > & v)
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<v.size(); ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result;

	make_array2d(result,n_rows,n_cols);

	for(size_t i=0;i<n_cols;++i)
	{
		for(size_t j=0;j<n_rows;++j)
		{
			result[j][i] = v[i][j];
		}
	}

	return result;
}

template<typename T>
std::vector< std::vector<T> > reverse_vertical(const std::vector< std::vector<T> > & v)
{
	size_t n_cols = v.size();
	size_t n_rows = 0;
	for(size_t i=0; i<v.size(); ++i)
		if(v[i].size() > n_rows) n_rows = v[i].size();

	std::vector< std::vector<T> > result;

	make_array2d(result,n_cols,n_rows);

	for(size_t i=0;i<n_cols;++i)
	{
		for(size_t j=0;j<n_rows;++j)
		{
			result[i][n_rows-j-1] = v[i][j];
		}
	}

	return result;
}

} // namespace brgastro

#endif // _BRG_MANIPULATIONS_HPP_INCLUDED_
