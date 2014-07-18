/**       @file SALTSA_vector.hpp
 *
 *     Project: SALTSA_example
 *  Repository: /SALTSA_example/SALTSA_vector.hpp
 *
 *  Created on: 24 Jun 2014
 *      Author: brg
 *
 *  This file contains element-wise operations for vectors, plus mathematical operations for vectors.
 */

#ifndef __SALTSA_VECTOR_HPP_INCLUDED__
#define __SALTSA_VECTOR_HPP_INCLUDED__

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "SALTSA_global.h"

namespace SALTSA {

// Element-wise functions
#if (1)

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
const std::vector<T1> apply( const f func, const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), func);

	return result;
}

template< typename f, typename T1, typename T2 >
const std::vector<T1> apply( const f func, const T1 & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++) result = func(v1,v2[i]);

	return result;
}

template< typename f, typename T1, typename T2 >
const std::vector<T1> apply( const f func, const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++) result = func(v1,v2[i]);

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
const std::vector<T> rand_vector_of_size(const f func, const int size)
{
	std::vector<T> result(size);

	for(unsigned int i = 0; i < size; i++) result[i] = func();

	return result;
}

template<typename f, typename T1>
const std::vector<T1> rand_vector_of_size(const f func, const T1 & v1, const int size)
{
	std::vector<T1> result(size);

	for(unsigned int i = 0; i < size; i++) result[i] = func(v1);

	return result;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector_of_size(const f func, const T1 & v1, const T2 & v2, const int size)
{
	std::vector<T1> result(size);

	for(unsigned int i = 0; i < size; i++) result[i] = func(v1,v2);

	return result;
}

template<typename f, typename T1>
const std::vector<T1> rand_vector(const f func, const std::vector<T1> & v1)
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++) result[i] = func(v1[i]);

	return result;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, const std::vector<T1> & v1, const std::vector<T2> & v2)
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++) result[i] = (func)(v1[i],v2.at(i));

	return result;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, const std::vector<T1> & v1, const T2 & v2)
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++) result[i] = func(v1[i],v2);

	return result;
}

template<typename f, typename T1, typename T2>
const std::vector<T1> rand_vector(const f func, const T1 & v1, const std::vector<T2> & v2)
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++) result[i] = func(v1,v2[i]);

	return result;
}

#endif

// Element-wise addition
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> add( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] + v2.at(i);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> add( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] + v2;
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> add( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = v2[i] + v1;
	}

	return result;
}

template< typename T1, typename T2 >
const T1 add( const T1 & v1, const T2 & v2 )
{
	return v1 + v2;
}

#endif // Element-wise addition

// Element-wise subtraction
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> subtract( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] - v2.at(i);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> subtract( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] - v2;
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> subtract( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = v1 - v2[i];
	}

	return result;
}

template< typename T1, typename T2 >
const T1 subtract( const T1 & v1, const T2 & v2 )
{
	return v1 - v2;
}

#endif // Element-wise subtraction

// Element-wise multiplication
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> multiply( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] * v2.at(i);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> multiply( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] * v2;
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> multiply( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = v2[i] * v1;
	}

	return result;
}

template< typename T1, typename T2 >
const T1 multiply( const T1 & v1, const T2 & v2 )
{
	return v1 * v2;
}

#endif // Element-wise multiplication

// Element-wise division
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> divide( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] / v2.at(i);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> divide( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = v1[i] / v2;
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> divide( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = v1 / v2[i];
	}

	return result;
}

template< typename T1, typename T2 >
const T1 divide( const T1 & v1, const T2 & v2 )
{
	return v1 / v2;
}

#endif // Element-wise divide

// Element-wise power
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> pow( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = std::pow(v1[i], v2.at(i));
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> pow( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = pow(v1[i], v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> pow( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = pow(v1, v2[i]);
	}

	return result;
}

template< typename T1, typename T2 >
const T1 pow( const T1 & v1, const T2 & v2 )
{
	return std::pow(v1, v2);
}

#endif // Element-wise power

// Element-wise max
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> max( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = max(v1[i], v2.at(i));
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> max( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = max(v1[i], v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> max( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = max(v1, v2[i]);
	}

	return result;
}

#endif // Element-wise max

// Element-wise min
#if (1)

template< typename T1, typename T2 >
const std::vector<T1> min( const std::vector<T1> & v1, const std::vector<T2> &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = min(v1[i], v2.at(i));
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> min( const std::vector<T1> & v1, const T2 &v2 )
{
	std::vector<T1> result(v1.size());

	for(unsigned int i = 0; i < v1.size(); i++)
	{
		result[i] = min(v1[i], v2);
	}

	return result;
}

template< typename T1, typename T2 >
const std::vector<T1> min( const T2 & v1, const std::vector<T1> &v2 )
{
	std::vector<T1> result(v2.size());

	for(unsigned int i = 0; i < v2.size(); i++)
	{
		result[i] = min(v1, v2[i]);
	}

	return result;
}

#endif // Element-wise min

// Element-wise negate
#if (1)

template< typename T >
const std::vector<T> negate( const std::vector<T> & v )
{
	std::vector<T> result(v.size(),0);

	for(unsigned int i = 0; i < v.size(); i++)
	{
		result[i] = -v[i];
	}

	return result;
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
const std::vector<T> abs( const std::vector<T> & v )
{
	std::vector<T> result(v.size(),0);

	for(unsigned int i = 0; i < v.size(); i++)
	{
		result[i] = std::abs(v[i]);
	}

	return result;
}

template< typename T >
const T abs( const T & v )
{
	return std::abs(v);
}

#endif // Element-wise abs

// Element-wise square root
#if (1)

template< typename T >
const std::vector<T> sqrt( const std::vector<T> & v )
{
	std::vector<T> result(v.size(),0);

	for(unsigned int i = 0; i < v.size(); i++)
	{
		result[i] = sqrt(v[i]);
	}

	return result;
}

template< typename T >
const T sqrt( const T & v )
{
	return std::sqrt(v);
}

#endif // Element-wise square root

// Element-wise exponential
#if (1)

template< typename T >
const std::vector<T> exp( const std::vector<T> & v )
{
	std::vector<T> result(v.size());

	for(unsigned int i = 0; i < v.size(); i++)
		result[i] = std::exp(v[i]);

	return result;
}

template< typename T >
const T exp( const T & v )
{
	return std::exp(v);
}

#endif // Element-wise exponential

// Element-wise safe_d
#if (1)

template< typename T >
const std::vector<T> safe_d( const std::vector<T> & v )
{
	std::vector<T> result(v.size());

	for(unsigned int i = 0; i < v.size(); i++)
		result[i] = safe_d(v[i]);

	return result;
}

#endif // Element-wise safe_d

#endif // Element-wise functions

// Mathematical functions
#if (1)

// Sum
#if (1)

template< typename T >
const T sum(const std::vector<T> v)
{
	T result = 0;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		result += v[i];
	}
	return result;
}

template< typename T >
const T sum(const T v)
{
	return v;
}

#endif // Sum

// Product
#if (1)

template< typename T >
const T product(const std::vector<T> v)
{
	T result = 1;
	for(unsigned int i = 0; i < v.size(); i++)
	{
		result *= v[i];
	}
	return result;
}

template< typename T >
const T product(const T v)
{
	return v;
}

#endif // Product

// Mean
#if (1)

template< typename T >
const T mean(const std::vector<T> v)
{
	if(v.size()==0) return 0;
	return sum(v)/v.size();
}

template< typename T >
const T mean(const T v)
{
	return v;
}

#endif // Mean

// Standard Deviation
#if (1)

template< typename T >
const T std(const std::vector<T> v)
{
	if(v.size()<=1) return 0;

	return sqrt( divide(subtract(sum( pow(v,2) ), pow(sum(v),2) ), v.size() ) );
}

template< typename T >
const T std(const T v)
{
	return 0;
}

#endif // Standard Deviation

#endif // Mathematical functions


} // namespace SALTSA




#endif // __SALTSA_VECTOR_HPP_INCLUDED__
