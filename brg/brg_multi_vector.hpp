/**       @file brg_multi_vector.hpp
 *
 *     Project: brg
 *        Path: /brg/brg_multi_vector.hpp
 *
 *  Created on: 16 Jul 2014
 *      Author: brg
 */

#ifndef __BRG_MULTI_VECTOR_HPP_INCLUDED__
#define __BRG_MULTI_VECTOR_HPP_INCLUDED__

#include <vector>

namespace brgastro {

template < typename T, unsigned int NumDims >
class multi_vector{
private:
	std::vector< multi_vector<T,NumDims-1> > base;

	T & _at(std::vector<unsigned int> & position) const
	{
		unsigned int i = position.back();
		if(i >= base.size())
			throw std::out_of_range;
		position.pop_back();
		return base[i]._at(position);
	}

	T & _(std::vector<unsigned int> & position) const
	{
		unsigned int i = position.back();
		position.pop_back();
		return base[i]._(position);
	}

	void _resize(std::vector<unsigned int> & new_size)
	{
		unsigned int i = new_size.back();
		new_size.pop_back();
		base.resize(i);
		for(int j = 0; j < i; j++)
		{
			base[j]._resize(new_size);
		}
	}

	void _resize(std::vector<unsigned int> & new_size, T new_value)
	{
		unsigned int i = new_size.back();
		new_size.pop_back();
		base.resize(i);
		for(int j = 0; j < i; j++)
		{
			base[j]._resize(new_size,new_value);
		}
	}

public:
	T & at(const std::vector<unsigned int> & position) const
	{
		std::vector<unsigned int> new_position = position;
		unsigned int i = new_position.back();
		if(i >= base.size())
			throw std::out_of_range;
		new_position.pop_back();
		return base[i]._at(new_position);
	}

	T & operator()(const std::vector<unsigned int> & position) const
	{
		std::vector<unsigned int> new_position = position;
		unsigned int i = new_position.back();
		new_position.pop_back();
		return base[i]._(new_position);
	}

	const T & operator[](unsigned int i) const
	{
		return base[i];
	}

	void resize(const std::vector<unsigned int> & new_size)
	{
		std::vector<unsigned int> new_new_size = new_size;
		unsigned int i = new_new_size.back();
		new_new_size.pop_back();
		base.resize(i);
		for(int j = 0; j < i; j++)
		{
			base[j]._resize(new_new_size);
		}
	}

	void resize(const std::vector<unsigned int> & new_size, T new_value)
	{
		std::vector<unsigned int> new_new_size = new_size;
		unsigned int i = new_new_size.back();
		new_new_size.pop_back();
		base.resize(i);
		for(int j = 0; j < i; j++)
		{
			base[j]._resize(new_new_size,new_value);
		}
	}

	void clear()
	{
		base.clear();
	}

	const std::vector<unsigned int> & shape() const
	{
		std::vector<unsigned int> res_shape = base[0].shape();
		res_shape.push_back(base.size());
		return res_shape;
	}
};

template < typename T>
class multi_vector<T,1>{
private:
	std::vector<T> base;
public:
	T & at(const std::vector<unsigned int> & position) const
	{
		return base.at(position.back());
	}

	T & operator()(const std::vector<unsigned int> & position) const
	{
		return base[position.back()];
	}

	T & operator[](unsigned int i) const
	{
		return base[i];
	}

	void resize(const std::vector<unsigned int> & new_size)
	{
		base.resize(new_size.back());
	}

	void resize(const std::vector<unsigned int> & new_size, T new_value)
	{
		base.resize(new_size, new_value);
	}

	void clear()
	{
		base.clear();
	}

	const std::vector<unsigned int> & shape() const
	{
		return std::vector<unsigned int>(1,base.size());
	}
};

template < typename T>
class multi_vector<T,0>{
private:
	T base;
public:
	T & at(const std::vector<unsigned int> & position)
	{
		return base;
	}

	T & operator()(const std::vector<unsigned int> & position)
	{
		return base;
	}

	const T & operator[](unsigned int i) const
	{
		return base;
	}

	void resize(const std::vector<unsigned int> & new_size)
	{
		throw std::runtime_error("ERROR: D=0 multi_vector cannot be resized.\n");
	}

	void resize(const std::vector<unsigned int> & new_size, T new_value)
	{
		throw std::runtime_error("ERROR: D=0 multi_vector cannot be resized.\n");
	}

	void clear()
	{
		set_zero(base);
	}

	const std::vector<unsigned int> & shape() const
	{
		return std::vector<unsigned int>();
	}
};



} // namespace brg_astro



#endif // __BRG_MULTI_VECTOR_HPP_INCLUDED__
