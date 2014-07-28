/**       @file brg_vector.hpp
 *
 *     Project: brg
 *        Path: /brg/brg_vector.hpp
 *
 *  Created on: 23 Jul 2014
 *      Author: brg
 */

#ifndef __BRG_VECTOR_HPP_INCLUDED__
#define __BRG_VECTOR_HPP_INCLUDED__

#include <vector>
#include <iterator>

#include "brg_vector_functions.hpp"

namespace brgastro {

template<typename T, typename A=T>
class vector: private std::vector<T,A>
{
	typedef std::vector<T,A> v;
	typedef std::vector<T,A>::size_type vsize_t;
	typedef unsigned short int dsize_t;
	typedef unsigned int ssize_t;
	typedef std::vector< unsigned int > shape_t;
private:
	dsize_t _num_dim_;
	shape_t _shape_;

	vsize_t _get_p(const shape_t & position) const
	{
		vsize_t p = 0;
		vsize_t m = 1;

		for(ssize_t i=0; i<position.size(); i++ )
		{
			p += m * position[i];
			m *= _shape_[i];
		}
		return p;
	}
public:

    // Using all vector functions except operator=, resize, assign, insert, erase, clear,
	// push_back, and pop_back
#if (1)
	using v::begin;
	using v::end;
	using v::rbegin;
	using v::rend;

	using v::size;
	using v::max_size;
	using v::capacity;
	using v::empty;
	using v::reserve;

	using v::operator[];
	using v::at;
	using v::front;
	using v::back;

	using v::assign;
	using v::push_back;
	using v::pop_back;
	using v::insert;
	using v::erase;
	using v::clear;

	using v::get_allocator;

#ifdef _BRG_USE_CPP_11_STD_
	using v::cbegin;
	using v::cend;
	using v::crbegin;
	using v::crend;

	using v::shrink_to_fit;

	using v::data;
#endif

#endif

	// Specialized versions of vector functions
#if (1)
	template <class InputIterator>
	void assign (InputIterator first, InputIterator last)
	{
		v::assign(first,last);
		reshape(shape_t(_num_dim_,size()));
	}
	void assign (v::size_type n, const T & val)
	{
		v::assign(n,val);
		reshape(shape_t(_num_dim_,n));
	}
	vector<T,A &> push_back (const T & val)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		_shape_[0] += 1;
		try
		{
			v::push_back(val);
		}
		catch(const std::exception &e)
		{
			_shape_[0] -= 1;
			throw;
		}
		return *this;
	}
	iterator insert (iterator position, const T & val)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		_shape_[0] += 1;
		try
		{
			return v::insert(position,val);
		}
		catch(const std::exception &e)
		{
			_shape_[0] -= 1;
			throw;
		}
	}
	void insert (iterator position, size_type n, const value_type& val)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		try
		{
			v::insert(position,n,val);
			_shape_[0] += n;
		}
		catch(const std::exception &e)
		{
			throw;
		}
	}
	template <class InputIterator>
	void insert (iterator position, InputIterator first, InputIterator last)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		try
		{
			v::insert(position,n,val);
			_shape_[0] = size();
		}
		catch(const std::exception &e)
		{
			throw;
		}
	}
	iterator erase (iterator position)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		_shape_[0] -= 1;
		try
		{
			return v::erase(position);
		}
		catch(const std::exception &e)
		{
			_shape_[0] += 1;
			throw;
		}
	}
	iterator erase (iterator first, iterator last)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		try
		{
			iterator result = v::erase(position);
			_shape_[0] = size();
			return result;
		}
		catch(const std::exception &e)
		{
			throw;
		}
	}
	void resize (vsize_t n, T val = T())
	{
		_num_dim_ = 1;
		_shape_(_num_dim_,n);
		v::resize(n,val);
	}
	const T & pop_back()
	{
		T result = v::back();
		v::pop_back();
		return T;
	}
	void del_back()
	{
		v::pop_back();
	}
	void clear()
	{
		v::clear();
		reshape(shape_t(0,0));
	}

#ifdef _BRG_USE_CPP_11_STD_
	template <class... Args>
	iterator emplace (const_iterator position, Args&&... args)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		_shape_[0] += 1;
		try
		{
			v::emplace(position,args);
		}
		catch(const std::exception &e)
		{
			_shape_[0] -= 1;
			throw;
		}
	}
	template <class... Args>
	iterator emplace_back (Args&&... args)
	{
		if(_num_dim_ != 1)
			throw std::out_of_range;
		_shape_[0] += 1;
		try
		{
			v::emplace_back(args);
		}
		catch(const std::exception &e)
		{
			_shape_[0] -= 1;
			throw;
		}
	}
#endif

#endif // Specialized versions of vector functions

	// Overloaded assignment(-related) operators
#if (1)
	vector<T, A> & operator=(vector<T, A> other)
	{
		this->swap(other);
		return *this;
	}

	vector<T, A> & copy(vector<T, A> other)
	{
		return *this = other;
	}

	vector<T, A> & swap(vector<T, A> & other)
	{
		v::swap(other);
		std::swap(_num_dim_,other._num_dim_);
		std::swap(_shape_,other._shape_);
		return *this;
	}

	operator T() const
	{
		return v::at(0);
	}


	vector<T, A> & operator+=( vector<T, A> other )
	{
		if(_shape_ != other._shape_)
			throw std::out_of_range;
		for(vsize_t i=0; i<size(); i++)
			v[i] += other[i];
		return *this;
	}
	template<typename T_o>
	vector<T, A> & operator+=( T_o other )
	{
		// Special handling for scalars
		for(vsize_t i=0; i<size(); i++)
			v[i] += other;
		return *this;
	}
	vector<T, A> & operator-=( vector<T, A> other )
	{
		if(_shape_ != other._shape_)
			throw std::out_of_range;
		for(vsize_t i=0; i<size(); i++)
			v[i] -= other[i];
		return *this;
	}
	template<typename T_o>
	vector<T, A> & operator-=( T_o other )
	{
		// Special handling for scalars
		for(vsize_t i=0; i<size(); i++)
			v[i] -= other;
		return *this;
	}
	vector<T, A> & operator*=( vector<T, A> other )
	{
		if(_shape_ != other._shape_)
			throw std::out_of_range;
		for(vsize_t i=0; i<size(); i++)
			v[i] *= other[i];
		return *this;
	}
	template<typename T_o>
	vector<T, A> & operator*=( T_o other )
	{
		// Special handling for scalars
		for(vsize_t i=0; i<size(); i++)
			v[i] *= other;
		return *this;
	}
	vector<T, A> & operator/=( vector<T, A> other )
	{
		if(_shape_ != other._shape_)
			throw std::out_of_range;
		for(vsize_t i=0; i<size(); i++)
			v[i] += other[i];
		return *this;
	}
	template<typename T_o>
	vector<T, A> & operator/=( T_o other )
	{
		// Special handling for scalars
		for(vsize_t i=0; i<size(); i++)
			v[i] += other;
		return *this;
	}

	vector<T, A> & operator++()
	{
		for(vsize_t i=0; i<size(); i++)
			v[i] += 1;
		return *this;
	}
	const vector<T, A> & operator++(int)
	{
		vector<T,A> old(*this);
		for(vsize_t i=0; i<size(); i++)
			v[i] += 1;
		return old;
	}

#endif

	// New methods and accessors
#if (1)

	// Accessors
#if (1)
	const shape_t & shape() const
	{
		return _shape_;
	}

	const dsize_t & num_dim() const
	{
		return _num_dim_;
	}

	// Unsafe element-access by position vector
	T & operator() (shape_t position)
	{
		return v[_get_p(position)];
	}
	const T & operator() (shape_t position) const
	{
		return v[_get_p(position)];
	}

	// Safe element-access by position vector
	T & at(shape_t position)
	{
		if(position.size() != _num_dim_)
			throw std::out_of_range;
		if(not_all_true(position<_shape_))
			throw std::out_of_range;

		return (*this)(position);
	}
	const T & at(shape_t position) const
	{
		if(position.size() != _num_dim_)
			throw std::out_of_range;
		if(not_all_true(position<_shape_))
			throw std::out_of_range;

		return (*this)(position);
	}
#endif

	// Reshaping
	void reshape(const shape_t & new_shape, const T & new_val = T())
	{
		// Check if no reshaping is necessary
		if(_shape_==new_shape) return;

		_shape_ = new_shape;
		if(_shape_.size() > std::numeric_limits<dsize_t>::max())
			throw std::out_of_range;
		_num_dim_ = _shape_.size();

		vsize_t new_size = 1;
		for(shape_t::iterator itr=new_shape.begin(); itr!=new_shape.end();++itr)
			new_size *= *itr;

		if(new_size >= std::vector<T,A>::max_size())
			throw std::out_of_range; // There's a slim chance this will throw on a valid size, but unlikely

		v::resize(new_size,new_val);
	}

#endif // New methods

	// Conversion to a std::vector
	operator std::vector<T,A>() const
	{
		return v;
	}

	// Constructors and destructors
#if (1)
	vector()
	: _num_dim_(0)
	{
	}
	vector( vsize_t init_size, T init_val )
	{
		resize(init_size,init_val);
	}
	vector( shape_t init_shape, T init_val )
	{
		reshape(init_shape,init_val);
	}
	template<typename T_o, typename A_o>
	vector(const vector<T_o, A_o> & other)
	{
		reshape(other._shape_);
		for(vsize_t i=0; i<other_vector.size(); i++) v[i] = other[i];
	}
	template<typename T_o, typename A_o>
	vector(const std::vector<T_o, A_o> & other)
	{
		resize(other.size());
		for(vsize_t i=0; i<other_vector.size(); i++) v[i] = other[i];
	}
	template<typename T_o>
	vector( const T (&array)[N] )
	{
		resize(N);
		for(vsize_t i=0; i<N; i++) v[i] = array[i];
	}
	template<typename T_o, typename A_o>
	vector(const T_o & other)
	{
		resize((vsize_t)1);
		v[0] = other;
	}
	virtual ~vector()
	{
		// Everything will be handled automatically by _shape_ and the data's destructors
	}

	// Move Constructors
#ifdef _BRG_USE_CPP_11_
	template<typename T_o, typename A_o>
	vector(std::vector<T_o, A_o> && other)
	: vector()
	{
		swap(*this,other);
	}
	template<typename T_o, typename A_o>
	vector(vector<T_o, A_o> && other)
	: vector()
	{
		swap(*this,other);
	}
#endif // _BRG_USE_CPP_11_

#endif

}; // class vector

// Global functions
#if (1)

template<typename T, typename A>
void swap(vector<T, A> & v1, vector<T, A> & v2)
{
	v1.swap(v2);
}

#endif // Global functions

// Overloaded operators
#if (1)

// Arithmetic operators
#if (1)

template<typename T1, typename A1, typename T2, typename A2>
const vector<T1, A1> & operator+( vector<T1, A1> v1, const vector<T2, A2> & v2 ) const
{
	v1 += v2;
	return v1;
}
template<typename T1, typename A1, typename T2>
const vector<T1, A1> & operator+( vector<T1, A1> v1, T2 v2 ) const
{
	v1 += v2;
	return v1;
}
template<typename T1, typename T2, typename A2>
const vector<T2, A2> & operator+( T1 v1, vector<T2, A2> v2 ) const
{
	v2 += v1;
	return v2;
}
template<typename T1, typename A1, typename T2, typename A2>
const vector<T1, A1> & operator-( vector<T1, A1> v1, const vector<T2, A2> & v2 ) const
{
	v1 -= v2;
	return v1;
}
template<typename T1, typename A1, typename T2>
const vector<T1, A1> & operator-( vector<T1, A1> v1, T2 v2 ) const
{
	v1 -= v2;
	return v1;
}
template<typename T1, typename T2, typename A2>
const vector<T2, A2> & operator-( T1 v1, vector<T2, A2> v2 ) const
{
	for(vector<T2, A2>::vsize_t i=0; i<v2.size(); i++)
		v2[i] = v1-v2[i];
	return v2;
}
template<typename T1, typename A1, typename T2, typename A2>
const vector<T1, A1> & operator*( vector<T1, A1> v1, const vector<T2, A2> & v2 ) const
{
	v1 *= v2;
	return v1;
}
template<typename T1, typename A1, typename T2>
const vector<T1, A1> & operator*( vector<T1, A1> v1, T2 v2 ) const
{
	v1 *= v2;
	return v1;
}
template<typename T1, typename T2, typename A2>
const vector<T2, A2> & operator*( T1 v1, vector<T2, A2> v2 ) const
{
	v2 *= v1;
	return v2;
}
template<typename T1, typename A1, typename T2, typename A2>
const vector<T1, A1> & operator/( vector<T1, A1> v1, const vector<T2, A2> & v2 ) const
{
	v1 /= v2;
	return v1;
}
template<typename T1, typename A1, typename T2>
const vector<T1, A1> & operator/( vector<T1, A1> v1, T2 v2 ) const
{
	v1 /= v2;
	return v1;
}
template<typename T1, typename T2, typename A2>
const vector<T2, A2> & operator/( T1 v1, vector<T2, A2> v2 ) const
{
	for(vector<T2, A2>::vsize_t i=0; i<v2.size(); i++)
		v2[i] = v1/v2[i];
	return v2;
}

#endif // Arithmetic operators

// Comparison operators
#if (1)

// Element-wise equal
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator==( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	if(v1.size()==1)
	{
		std::vector<bool> result(v2.size());
		// Special handling for scalars
		T1 val = v1[0];
		for(brgastro::vector::vsize_t i=0; i<v2.size(); i++)
			result[i] = (val == v2[i]);
		return result;
	}
	else if(v2.size()==1)
	{
		std::vector<bool> result(v1.size());
		// Special handling for scalars
		T1 val = v2[0];
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			result[i] = (v1[i] == val);
		return result;
	}
	else
	{
		if(v1.shape()!= v2.shape())
			throw std::out_of_range;
		std::vector<bool> result(v1.size());
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			result[i] = (v1[i] == v2[i]);
		return result;
	}
}

// Element-wise not equal
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator!=( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	return !operator==(v1,v2);
}

// Element-wise less-than
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator<( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	if(v1.size()==1)
	{
		std::vector<bool> result(v2.size());
		// Special handling for scalars
		T1 val = v1[0];
		for(brgastro::vector::vsize_t i=0; i<v2.size(); i++)
			result[i] = (val < v2[i]);
		return result;
	}
	else if(v2.size()==1)
	{
		std::vector<bool> result(v1.size());
		// Special handling for scalars
		T1 val = v2[0];
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			result[i] = (v1[i] < val);
		return result;
	}
	else
	{
		if(v1.shape()!= v2.shape())
			throw std::out_of_range;
		std::vector<bool> result(v1.size());
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			result[i] = (v1[i] < v2[i]);
		return result;
	}
}

// Element-wise greater-than
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator>( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	return operator<(v2,v1);
}

// Element-wise less-than or equal to
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator<=( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	return !operator>(v1,v2);
}

// Element-wise greater-than or equal to
template<typename T1, typename A1, typename T2, typename A2>
const std::vector<bool> & operator>=( const vector<T1, A1> & v1, const vector<T2, A2> & v2 ) const
{
	return !operator<(v1,v2);
}

#endif // Comparison operators

#endif // Overloaded operators

// Element-wise min/max/bound
#if (1)

template<typename T, typename A, typename T_o, typename A_o>
const vector<T,A> & min( vector<T,A> v1, const vector<T_o,A_o> & v2 )
{
	if(v2.size()==1)
	{
		// Special handling for scalars
		T val = v2[0];
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			v1[i] = min(v1[i],val);
		return v1;
	}
	else if(v1.size()==1)
	{
		// Special handling for scalars
		const vector<T,A> result(v2);
		T val = v1[0];
		for(brgastro::vector::vsize_t i=0; i<v2.size(); i++)
			result[i] = min(val,result[i]);
		return result;
	}
	else
	{
		if(v1.shape()!= v2.shape())
			throw std::out_of_range;
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			v1[i] = min(v1[i],v2[i]);
		return v1;
	}
}

template<typename T, typename A, typename T_o, typename A_o>
const vector<T,A> & max( vector<T,A> v1, const vector<T_o,A_o> & v2 )
{
	if(v2.size()==1)
	{
		// Special handling for scalars
		T val = v2[0];
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			v1[i] = max(v1[i],val);
		return v1;
	}
	else if(v1.size()==1)
	{
		// Special handling for scalars
		const vector<T,A> result(v2);
		T val = v1[0];
		for(brgastro::vector::vsize_t i=0; i<v2.size(); i++)
			result[i] = max(val,result[i]);
		return result;
	}
	else
	{
		if(v1.shape()!= v2.shape())
			throw std::out_of_range;
		for(brgastro::vector::vsize_t i=0; i<v1.size(); i++)
			v1[i] = max(v1[i],v2[i]);
		return v1;
	}
}

template<typename T_1, typename A_1, typename T_2, typename A_2, typename T_3, typename A_3>
const vector<T_2,A_2> & bound( const vector<T_1,A_1> & v1, const vector<T_2,A_2> & v2, const vector<T_3,A_3> & v3 )
{
	return min(max(v2,v1),v3);
}

#endif // Element-wise min/max/bound

// Math and safe functions
#if (1)

// Element-wise power
#if (1)

template< typename T1, typename T2 >
const vector<T1> pow( vector<T1> v1, const vector<T2> &v2 )
{
	if(v1.shape()!= v2.shape())
		throw std::out_of_range;
	for(unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] = std::pow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const vector<T1> pow( const vector<T1> & v1, const T2 &v2 )
{
	for(unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] = std::pow(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const vector<T1> pow( const T2 & v1, vector<T1> v2 )
{
	for(unsigned int i = 0; i < v2.size(); i++)
	{
		v2[i] = std::pow(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise power

// Element-wise safe power
#if (1)

template< typename T1, typename T2 >
const vector<T1> safe_pow( vector<T1> v1, const vector<T2> &v2 )
{
	if(v1.shape()!= v2.shape())
		throw std::out_of_range;
	for(unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] = safe_pow(v1[i], v2[i]);
	}

	return v1;
}

template< typename T1, typename T2 >
const vector<T1> safe_pow( const vector<T1> & v1, const T2 &v2 )
{
	for(unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] = safe_pow(v1[i], v2);
	}

	return v1;
}

template< typename T1, typename T2 >
const vector<T1> safe_pow( const T2 & v1, vector<T1> v2 )
{
	for(unsigned int i = 0; i < v2.size(); i++)
	{
		v2[i] = safe_pow(v1, v2[i]);
	}

	return v2;
}

#endif // Element-wise safe power

// Element-wise negate
#if (1)

template< typename T >
const vector<T> operator-( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = -v[i];
	}
	return v;
}

#endif // Element-wise negate

// Element-wise abs
#if (1)

template< typename T >
const vector<T> abs( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = std::abs(v[i]);
	}

	return v;
}

#endif // Element-wise abs

// Element-wise square root
#if (1)

template< typename T >
const vector<T> sqrt( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = std::sqrt(v[i]);
	}

	return v;
}

#endif // Element-wise square root

// Element-wise safe square root
#if (1)

template< typename T >
const vector<T> safe_sqrt( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = safe_sqrt(v[i]);
	}

	return v;
}

#endif // Element-wise safe square root

// Element-wise exponential
#if (1)

template< typename T >
const vector<T> exp( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
		v[i] = std::exp(v[i]);

	return v;
}

#endif // Element-wise exponential

// Element-wise safe_d
#if (1)

template< typename T >
const vector<T> safe_d( vector<T> v )
{
	for(unsigned int i = 0; i < v.size(); i++)
		v[i] = safe_d(v[i]);

	return v;
}

#endif // Element-wise safe_d

#endif // Math and safe functions

} // namespace brgastro

#endif // __BRG_VECTOR_HPP_INCLUDED__
