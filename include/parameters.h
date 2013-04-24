#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <stdexcept>
#include <cstddef>

class BadParameters : public std::runtime_error
{
public:
	BadParameters(const std::string& what) : std::runtime_error("Parameters error: " + what) {}
};

/// A fifo buffer for arbitrary parameter data.
class Parameters
{
public:
	Parameters()
	: cur(0)
	{
	}
	
	Parameters(const Parameters& other)
	{
		values.reserve(other.values.size());
		// clone all contained values
		for(std::size_t i = 0, n = other.values.size(); i < n; ++i)
			values.push_back(other.values[i]->clone());
		cur = other.cur;
	}
	
	~Parameters()
	{
		this->clear();
	}
	
	Parameters& operator=(Parameters rhs)
	{
		swap(*this, rhs);
		return *this;
	}
	
	/// Delete all contained values
	void clear()
	{
		for(std::size_t i = 0, n = values.size(); i < n; ++i)
			delete values[i];
		values.clear();
		cur = 0;
	}
	
	/// Check if parameters are present.
	bool empty()
	{
		return values.empty();
	}
	
	/// Reset to beginning of parameters list.
	void reset()
	{
		cur = 0;
	}
	
	/// Input
	template<typename T>
	friend Parameters& operator<<(Parameters& p, const T& x)
	{
		p.values.push_back(new value<T>(x));
		return p;
	}
	
	/// Output
	template<typename T>
	friend Parameters& operator>>(Parameters& p, T& x)
	{
		if(p.cur == p.values.size())
			throw BadParameters("Reading from past the end.");
		
		value<T>* v = dynamic_cast<value<T>*>(p.values[p.cur++]);
		
		if(v == 0)
			throw BadParameters("Could not cast parameters to requested type.");
		
		x = *v;
		
		return p;
	}
	
	friend void swap(Parameters& a, Parameters& b)
	{
		using std::swap;
		swap(a.cur, b.cur);
		swap(a.values, b.values);
	}
	
private:
	/// base class for type hiding
	class base_value
	{
	public:
		virtual ~base_value() {}
		virtual base_value* clone() = 0;
	};
	
	/// value class
	template<typename value_type>
	class value : public base_value
	{
	public:
		value(const value_type& x) : v(x) {}
		virtual ~value() {}
		
		operator const value_type&() { return v; }
		
		value<value_type>* clone() { return new value<value_type>(v); }
		
	private:
		const value_type v;
	};
	
	std::vector<base_value*> values;
	std::size_t cur;
};

#endif
