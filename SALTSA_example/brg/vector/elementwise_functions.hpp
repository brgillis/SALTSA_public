/**********************************************************************\
  @file elementwise_functions.hpp

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

#ifndef _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_
#define _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "brg/global.h"

#include "brg/math/misc_math.hpp"
#include "brg/math/safe_math.hpp"

namespace brgastro {

// These apply a standard function to each element of a vector and return a vector of results

// Element-wise generic function
#if (1)

template< typename f, typename T >
const std::vector<T> apply( const f func, const std::vector<T> & v1 )
{
	std::vector<T> result(v1.size());

	std::transform(v1.begin(), v1.end(), result.begin(), func);

	return result;
}

template< typename f, typename T1, typename T2 >
const std::vector<T1> apply( const f func, const std::vector<T1> & v1, const std::vector<T2> & v2 )
{
	std::vector<T1> result(v1.size());

	std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), func);

	return result;
}

template< typename f, typename T1, typename T2 >
const std::vector<T1> apply( const f func, const T1 & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++) result = func(v1,v2[i]);

	return result;
}

template< typename f, typename T1, typename T2 >
const std::vector<T1> apply( const f func, const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++) result = func(v1,v2[i]);

	return result;
}

template< typename f, typename T1, typename T2 >
const T1 apply( const f * func, const T1 & v1, const T2 &v2 )
{
	T1 result = func(v1,v2);

	return result;
}

#endif // Element-wise generic function

// Random values
#if (1)

template<typename T, typename f>
const std::vector<T> rand_vector_of_size(const f func, const size_t size)
{
	std::vector<T> result(size);

	for(size_t i = 0; i < size; i++) result[i] = func();

	return result;
}

template<typename f, typename T1>
const std::vector<T1> rand_vector_of_size(const f func, const T1 & v1, const size_t size)
{
	std::vector<T1> result(size);

	for(size_t i = 0; i < size; i++) result[i] = func(v1);

	return result;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector_of_size(const f func, const T1 & v1, const T2 & v2, const size_t size)
{
	std::vector<T1> result(size);

	for(size_t i = 0; i < size; i++) result[i] = func(v1,v2);

	return result;
}

template<typename f, typename T1>
const std::vector<T1> rand_vector(const f func, std::vector<T1> v1)
{
	for(size_t i = 0; i < v1.size(); i++) v1[i] = func(v1[i]);

	return v1;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, std::vector<T1> v1, const std::vector<T2> & v2)
{
	assert(v1.size()==v2.size());

	for(size_t i = 0; i < v1.size(); i++) v1[i] = (func)(v1[i],v2[i]);

	return v1;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, std::vector<T1> v1, const T2 & v2)
{
	for(size_t i = 0; i < v1.size(); i++) v1[i] = func(v1[i],v2);

	return v1;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, const T1 & v1, const std::vector<T2> & v2)
{
	std::vector<T1> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++) result[i] = func(v1,v2[i]);

	return result;
}

#endif

// Element-wise math
#if (1)

// Element-wise addition
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> add( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] += v2[i];
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> add( std::vector<T1> v1, const T2 &v2 )
{
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] += v2;
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> add( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] += v1;
	}

	return v2;
}

template< typename T1, typename T2 >
const T1 add( T1 v1, const T2 & v2 )
{
	return v1 += v2;
}

#endif // Element-wise addition

// Element-wise subtraction
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> subtract( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] -= v2[i];
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> subtract( std::vector<T1> v1, const T2 &v2 )
{
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] -= v2;
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> subtract( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = v1-v2[i];
	}

	return v2;
}

template< typename T1, typename T2 >
const T1 subtract( T1 v1, const T2 & v2 )
{
	return v1 -= v2;
}

#endif // Element-wise subtraction

// Element-wise multiplication
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> multiply( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] *= v2[i];
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> multiply( std::vector<T1> v1, const T2 &v2 )
{
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] *= v2;
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> multiply( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] *= v1;
	}

	return v2;
}

template< typename T1, typename T2 >
const T1 multiply( T1 v1, const T2 & v2 )
{
	return v1 *= v2;
}

#endif // Element-wise multiplication

// Element-wise division
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> divide( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] /= v2[i];
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> divide( std::vector<T1> v1, const T2 &v2 )
{

	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] /= v2;
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> divide( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = v1/v2[i];
	}

	return v2;
}

template< typename T1, typename T2 >
const T1 divide( T1 v1, const T2 & v2 )
{
	return v1 /= v2;
}

#endif // Element-wise divide

// Element-wise power
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> pow( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = pow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> pow( std::vector<T1> v1, const T2 &v2 )
{
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = pow(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> pow( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = pow(v1, v2[i]);
	}

	return v2;
}

template< typename T1, typename T2 >
const T1 pow( const T1 & v1, const T2 & v2 )
{
	return std::pow(v1, v2);
}

template< typename T1, typename T2 >
const std::vector<T1> ipow( std::vector<T1> v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = ipow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> ipow( std::vector<T1> v1, const T2 &v2 )
{
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = ipow(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> ipow( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = ipow(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise power

// Element-wise safe power
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> safe_pow( std::vector<T1> v1, const std::vector<T2> &v2 )
{

	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = safe_pow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> safe_pow( std::vector<T1> v1, const T2 &v2 )
{

	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = safe_pow(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> safe_pow( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = safe_pow(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise safe power

// Element-wise max
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> max( std::vector<T1> v1, const std::vector<T2> &v2 )
{

	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = max(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> max( std::vector<T1> v1, const T2 &v2 )
{

	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = max(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> max( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = max(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise max

// Element-wise min
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> min( std::vector<T1> v1, const std::vector<T2> &v2 )
{

	assert(v1.size()==v2.size());
	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = min(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T1> min( std::vector<T1> v1, const T2 &v2 )
{

	for(size_t i = 0; i < v1.size(); i++)
	{
		v1[i] = min(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const std::vector<T2> min( const T1 & v1, std::vector<T2> v2 )
{
	for(size_t i = 0; i < v2.size(); i++)
	{
		v2[i] = min(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise min

// Element-wise negate
#if (1)

template< typename T >
const std::vector<T> negate( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = negate(v[i]);
	}

	return v;
}

template< typename T >
const T negate( const T & v )
{
	return -v;
}

#endif // Element-wise negate

// Element-wise abs
#if (1)

template< typename T >
const std::vector<T> abs( std::vector<T> v )
{
	using std::abs;
	using brgastro::abs;

	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = abs(v[i]);
	}

	return v;
}

#endif // Element-wise abs

// Element-wise square root
#if (1)

template< typename T >
const std::vector<T> sqrt( std::vector<T> v )
{
	using std::sqrt;
	using brgastro::sqrt;

	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = sqrt(v[i]);
	}

	return v;
}

#endif // Element-wise square root

// Element-wise safe square root
#if (1)

template< typename T >
const std::vector<T> safe_sqrt( std::vector<T> v )
{
	std::vector<T> result(v.size(),0);

	for(size_t i = 0; i < v.size(); i++)
	{
		result[i] = safe_sqrt(v[i]);
	}

	return result;
}

#endif // Element-wise safe square root

// Element-wise exponential
#if (1)

template< typename T >
const std::vector<T> exp( std::vector<T> v )
{
	using std::exp;
	using brgastro::exp;

	for(size_t i = 0; i < v.size(); i++)
		v[i] = exp(v[i]);

	return v;
}

#endif // Element-wise exponential

// Element-wise square
#if (1)

template< typename T >
const std::vector<T> square( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = square(v[i]);
	}

	return v;
}

#endif // Element-wise square

// Element-wise cube
#if (1)

template< typename T >
const std::vector<T> cube( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = cube(v[i]);
	}

	return v;
}

#endif // Element-wise cube

// Element-wise quart
#if (1)

template< typename T >
const std::vector<T> quart( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = quart(v[i]);
	}

	return v;
}

#endif // Element-wise quart

// Element-wise inverse
#if (1)

template< typename T >
const std::vector<T> inverse( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = inverse(v[i]);
	}

	return v;
}

#endif // Element-wise inverse

// Element-wise inv_square
#if (1)

template< typename T >
const std::vector<T> inv_square( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = inv_square(v[i]);
	}

	return v;
}

#endif // Element-wise inv_square

// Element-wise inv_cube
#if (1)

template< typename T >
const std::vector<T> inv_cube( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = inv_cube(v[i]);
	}

	return v;
}

#endif // Element-wise inv_cube

// Element-wise inv_quart
#if (1)

template< typename T >
const std::vector<T> inv_quart( std::vector<T> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = inv_quart(v[i]);
	}

	return v;
}

#endif // Element-wise inv_quart

// Element-wise safe_d
#if (1)

template< typename T >
const std::vector<T> safe_d( const std::vector<T> & v )
{
	std::vector<T> result(v.size());

	for(size_t i = 0; i < v.size(); i++)
		result[i] = safe_d(v[i]);

	return result;
}

#endif // Element-wise safe_d

#endif // Element-wise math

// Element-wise comparison (always return a vector of bools)
#if (1)

// Element-wise not
#if (1)

inline const std::vector<bool> operator!( std::vector<bool> v )
{
	for(size_t i = 0; i < v.size(); i++)
	{
		v[i] = !v[i];
	}
	return v;
}

#endif // Element-wise not

// Element-wise equal
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> equal( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] == v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> equal( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] == v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> equal( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v2[i] == v1);
	}

	return result;
}

template< typename T1, typename T2 >
const bool equal( const T1 & v1, const T2 & v2 )
{
	return (v1 == v2);
}

#endif // Element-wise equal

// Element-wise not_equal
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> not_equal( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] != v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> not_equal( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] != v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> not_equal( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v2[i] != v1);
	}

	return result;
}

template< typename T1, typename T2 >
const bool not_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 != v2);
}

#endif // Element-wise not_equal

// Element-wise less_than
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> less_than( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] < v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> less_than( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] < v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> less_than( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 < v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const bool less_than( const T1 & v1, const T2 & v2 )
{
	return (v1 < v2);
}

#endif // Element-wise less_than

// Element-wise greater_than
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> greater_than( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] > v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> greater_than( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] > v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> greater_than( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 > v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const bool greater_than( const T1 & v1, const T2 & v2 )
{
	return (v1 > v2);
}

#endif // Element-wise equal

// Element-wise less_than_or_equal
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> less_than_or_equal( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] <= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> less_than_or_equal( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] <= v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> less_than_or_equal( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 <= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const bool less_than_or_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 <= v2);
}

#endif // Element-wise less_than_or_equal

// Element-wise greater_than_or_equal
#if (1)
template< typename T1, typename T2 >
const std::vector<bool> greater_than_or_equal( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	assert(v1.size()==v2.size());

	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] >= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> greater_than_or_equal( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<bool> result(v1.size());

	for(size_t i = 0; i < v1.size(); i++)
	{
		result[i] = (v1[i] >= v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<bool> greater_than_or_equal( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<bool> result(v2.size());

	for(size_t i = 0; i < v2.size(); i++)
	{
		result[i] = (v1 >= v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const bool greater_than_or_equal( const T1 & v1, const T2 & v2 )
{
	return (v1 >= v2);
}

#endif // Element-wise greater_than_or_equal

#endif // Element-wise comparison

// Functions on bool vectors
#if (1)

// not
#if (1)
template<typename T>
inline const std::vector<T> v_not(const std::vector<T> v)
{
	for(size_t i=0; i < v.size(); i++)
	{
		v[i] = v_not(v[i]);
	}
	return v;
}

inline const bool v_not(const bool v)
{
	return !v;
}
#endif // not

#endif // Functions on bool vectors


} // namespace brgastro

#endif // _BRG_ELEMENTWISE_FUNCTIONS_HPP_INCLUDED_
