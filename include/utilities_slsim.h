#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "standard.h"
#include <typeinfo>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <iterator>
#include <cstdlib>
#include <random>
#if __cplusplus >= 201103L
#include <typeindex>
#endif

namespace Utilities
{
	// this is not for the user
	namespace detail
	{
#if __cplusplus < 201103L
		class type_index
		{
		public:
			type_index(const std::type_info& type) : t(type) {}
			inline bool operator<(const type_index& rhs) const { return t.before(rhs.t); }
			
		private:
			const std::type_info& t;
		};
#else
		using std::type_index;
#endif
	}
	
	template <class T>
	void Matrix(T **matrix, long rows, long cols){
		matrix = new T*[rows];
		matrix[0] = new T[rows*cols];
		for (long i = 1; i < rows; ++i)
			matrix[i] = matrix[0] + i * cols;
	}
	
	template <class T>
	void free_Matrix(T **matrix, long rows, long){
		if (rows) delete[] matrix[0];
		delete []matrix;
	}
	
	PosType **PosTypeMatrix(long rows, long cols);
	void free_PosTypeMatrix(PosType **matrix, long rows, long cols);
	
	/** \brief A container that can hold mixed objects all derived from
	 * a base class and retains the ability to access derived class functions/members.
	 *
	 * The indexing operator is overridden to return a reference to the base class and
	 * a templated index<>() function provides a reference to elements with derived class
	 * type information.
	 *
	 *
	 *<pre>
	 * example:
	 * InputParams params(paramfile);
     * params.put("gauss_r2",1,"non");
 	 * SourceGaussian source1(params);
	 * SourceUniform source2(params);
     *
	 * Utilities::MixedVector<Source> mvector;
     *
	 * mvector.push_back(source1);
	 * mvector.push_back(source1);
	 * mvector.push_back(source2);
	 * mvector.push_back(source2);
	 * std::cout << "Number of Uniform Sources " << mvector.size<SourceUniform>() << "   Number of Gausssian Sources "
	 *           << mvector.size<SourceGaussian>() << std::endl;
     *
     * // change derived class attribute
	 * mvector.get<SourceGaussian>(0).source_gauss_r2 = 0.5;
	 * std::cout << "A base class attribute " << mvector[2].getTotalFlux()
	 *           << " A derived class attribute  "
	 *           << mvector.get<SourceGaussian>(0).source_gauss_r2
	 *           << "  " << mvector.get<SourceGaussian>(1).source_gauss_r2
	 *           << std::endl;
	 * 
	 * // iterate all sources
	 * for(Utilities::MixedVector<Source>::iterator<> it = mvector.begin(); it != mvector.end(); ++it)
	 *   std::cout << "A base class attribute " << it->getTotalFlux() << std::endl;
	 * 
	 * // iterate SersicSources
	 * for(Utilities::MixedVector<Source>::iterator<SersicSource> it = mvector.begin<SersicSource>(); it != mvector.end<SersicSource>(); ++it)
	 *   std::cout << "A derived class attribute " << it->getSersicIndex() << std::endl;
	 *</pre>
	 */
	template<typename BaseT>
	class MixedVector
	{
	public: /* iterators */
		template<typename ValueT = BaseT> // TODO: needs to check subclass type
		class iterator
		{
		private:
			typedef typename std::vector<BaseT*>::iterator base_iterator;
			base_iterator it;
			
		public:
			typedef ValueT value_type;
			typedef ValueT* pointer;
			typedef ValueT& reference;
			typedef typename base_iterator::difference_type difference_type;
			typedef typename base_iterator::iterator_category iterator_category;
			
			iterator() {}
			iterator(base_iterator i) : it(i) {}
			iterator(const iterator& other) : it(other.it) {}
			
			iterator& operator=(const iterator& rhs) { it = rhs.it; return *this; }
			
			reference operator*() { return (reference)(**it); }
			const reference operator*() const { return (const reference)(**it); }
			
			pointer operator->() { return (pointer)(*it); }
			const pointer operator->() const { return (const pointer)(*it); }
			
			iterator& operator++() { ++it; return *this; }
			iterator operator++(int) { iterator tmp(*this); ++it; return tmp; }
			iterator& operator--() { --it; return *this; }
			iterator operator--(int) { iterator tmp(*this); --it; return tmp; }
			
			bool operator==(const iterator& rhs) const { return (it == rhs.it); }
			bool operator!=(const iterator& rhs) const { return (it != rhs.it); }
			
			iterator& operator+=(difference_type n) { it += n; return *this; }
			iterator& operator-=(difference_type n) { it -= n; return *this; }
			
			reference operator[](difference_type n) { return (reference)*it[n]; }
			const reference operator[](difference_type n) const { return (const reference)*it[n]; }
			
			friend iterator operator+(const iterator& i, difference_type n) { return iterator(i.it + n); }
			friend iterator operator+(difference_type n, const iterator& i) { return iterator(i.it + n); }
			friend iterator operator-(const iterator& i, difference_type n) { return iterator(i.it - n); }
			friend iterator operator-(difference_type n, const iterator& i) { return iterator(i.it - n); }
			
			friend difference_type operator-(const iterator& b, const iterator& a) { return (b.it - a.it); }
			
			friend bool operator<(const iterator&a, const iterator& b) { return (a.it < b.it); }
			friend bool operator>(const iterator&a, const iterator& b) { return (a.it > b.it); }
			friend bool operator<=(const iterator&a, const iterator& b) { return (a.it <= b.it); }
			friend bool operator>=(const iterator&a, const iterator& b) { return (a.it >= b.it); }
		};
		
	public:
		/// default constructor
		MixedVector()
		{
		}
		
		/// destroy the vector and all contained items
		~MixedVector()
		{
			clear();
		}
		
		/// swap two MixedVectors
		friend void swap(MixedVector<BaseT>& a, MixedVector<BaseT>& b)
		{
			using std::swap;
			
			swap(a.items, b.items);
			swap(a.tmap, b.tmap);
		}
		
		/// assign one MixedVector to another
		MixedVector<BaseT>& operator=(MixedVector<BaseT> rhs)
		{
			// copy-and-swap idiom
			swap(*this, rhs);
			return *this;
		}
		
		/// add an object of type SubclassT to the vector
		template<typename SubclassT>
		void push_back(const SubclassT& obj)
		{
			// make sure this is a subclass of BaseT
			check_subclass(obj);
			
			// copy the object
			SubclassT* copy = new SubclassT(obj);
			
			// add the copy of the object to the list of items
			items.push_back(copy);
			
			// add the copy to type map
			tmap[typeid(SubclassT)].push_back(copy);
		}
		
		/// pop element from back of vector
		void pop_back()
		{
			// get very last element
			BaseT* back = items.back();
			
			// remove from the vector in the type map
			tmap[back->type()].pop_back();
			
			// remove from items
			items.pop_back();
			
			// delete from memory
			delete back;
		}
		
		/// erase the last element of a specific type
		template<typename SubclassT>
		void pop_back()
		{
			// search for vector belonging to SubclassT
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return;
			
			// pop from vector for type
			found->second.pop_back();
			
			// go through list of all items from back, until item of correct type is found
			for(typename base_container::reverse_iterator it = items.rbegin(); it != items.rend(); ++it)
			{
				if((*it)->type() == typeid(SubclassT))
				{
					// remove the item
					delete (*it);
					items.erase(it);
					return;
				}
			}
		}
		
		/// clear all elements
		void clear()
		{
			for(std::size_t i = 0, n = items.size(); i < n; ++i)
				delete items[i];
			items.clear();
			tmap.clear();
		}
		
		/// clear all elements of a given type
		template<typename SubclassT>
		void clear()
		{
			// search for vector belonging to SubclassT
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return;
			
			base_container& titems = found.second;
						
			// delete items from vector for subclass
			for(std::size_t i = 0, n = titems.size(); i < n; ++i)
				delete titems[i];
			
			// remove items for subclass
			tmap.erase(found);
			
			// erase items matching type
			items.erase(std::remove_if(items.begin(), items.end(), is_type<SubclassT>), items.end());
		}
		
		/// direct access to underlying array
		BaseT* data()
		{
			return items.data();
		}
		
		/// direct access to underlying array (const)
		const BaseT* data() const
		{
			return items.data();
		}
		
		/// Checks if element i is of the derived type SubclassT
		template<typename SubclassT>
		bool type(std::size_t i)
		{
			return is_type<SubclassT>(items[i]);
		}
		
		/// Indexing operator for all elements.
		BaseT& operator[](std::size_t i)
		{
			return *items[i];
		}
		
		/// Indexing operator for all elements (const).
		const BaseT& operator[](std::size_t i) const
		{
			return *items[i];
		}
		
		/// indexed access
		BaseT& get(std::size_t i)
		{
			return *items[i];
		}
		
		/// indexed access (const)
		const BaseT& get(std::size_t i) const
		{
			return *items[i];
		}
		
		/// indexed access for type SubclassT
		template<typename SubclassT>
		SubclassT& get(std::size_t i)
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (SubclassT&)*found->second[i];
		}
		
		/// indexed access for type SubclassT (const)
		template<typename SubclassT>
		const SubclassT& get(std::size_t i) const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (const SubclassT&)*found->second[i];
		}
		
		/// indexed access with bounds checking
		BaseT& at(std::size_t i)
		{
			return *items.at(i);
		}
		
		/// indexed access with bounds checking (const)
		const BaseT& at(std::size_t i) const
		{
			return *items.at(i);
		}
		
		/** \brief Templated indexing operator for elements of a specific derived class.
		 * The index runs 0 ... size<SubclassT>() - 1
		 * Does bounds checking.
		 */
		template<typename SubclassT>
		SubclassT& at(std::size_t i)
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (SubclassT&)*found->second.at(i);
		}
		
		/// Templated indexing operator for elements of a specific derived class (const).
		template<typename SubclassT>
		const SubclassT& at(std::size_t i) const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (const SubclassT&)*found->second.at(i);
		}
		
		/// number of all elements
		std::size_t size() const
		{
			return items.size();
		}
		
		/// number of elements of a specific type
		template<typename SubclassT>
		std::size_t size() const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return 0;
			return found->second.size();
		}
		
		/// check if vector of items is empty
		bool empty() const
		{
			return items.empty();
		}
		
		/// check if vector of items of type SubclassT is empty
		template<typename SubclassT>
		bool empty() const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return true;
			
			return found->second.empty();
		}
		
		/// get iterator to first of all items
		iterator<> begin()
		{
			return iterator<>(items.begin());
		}
		
		/// get iterator to last of all items
		iterator<> end()
		{
			return iterator<>(items.end());
		}
		
		/// get iterator to first of items of type SubclassT
		template<typename SubclassT>
		iterator<SubclassT> begin()
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return iterator<SubclassT>(items.end());
			return iterator<SubclassT>(found->second.begin());
		}
		
		/// get iterator to last of items of type SubclassT
		template<typename SubclassT>
		iterator<SubclassT> end()
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return iterator<SubclassT>(items.end());
			return iterator<SubclassT>(found->second.end());
		}
		
	private:
		/// cannot copy MixedVector
		MixedVector(const MixedVector<BaseT>& other);
		
		typedef std::vector<BaseT*> base_container;
		
		inline static void check_subclass(const BaseT&) {}
		
		template<typename SubclassT>
		inline static bool is_type(const BaseT* obj) { return (typeid(*obj) == typeid(SubclassT)); }
		
		base_container items;
		
		typedef std::map<detail::type_index, base_container> tmap_t;
		tmap_t tmap;
	};
	
	/// A MixedVector for pointers.
	template<typename BaseT>
	class MixedVector<BaseT*>
	{
	public: /* iterators */
		template<typename ValueT = BaseT*> // TODO: needs to check subclass pointer type
		class iterator
		{
		private:
			typedef typename std::vector<BaseT*>::iterator base_iterator;
			base_iterator it;
			
		public:
			typedef typename std::vector<ValueT>::iterator::value_type value_type;
			typedef typename std::vector<ValueT>::iterator::pointer pointer;
			typedef typename std::vector<ValueT>::iterator::reference reference;
			typedef typename base_iterator::difference_type difference_type;
			typedef typename base_iterator::iterator_category iterator_category;
			
			iterator() {}
			iterator(base_iterator i) : it(i) {}
			iterator(const iterator& other) : it(other.it) {}
			
			iterator& operator=(const iterator& rhs) { it = rhs.it; return *this; }
			
			reference operator*() { return (reference)*it; }
			const reference operator*() const { return (const reference)*it; }
			
			iterator& operator++() { ++it; return *this; }
			iterator operator++(int) { iterator tmp(*this); ++it; return tmp; }
			iterator& operator--() { --it; return *this; }
			iterator operator--(int) { iterator tmp(*this); --it; return tmp; }
			
			bool operator==(const iterator& rhs) const { return (it == rhs.it); }
			bool operator!=(const iterator& rhs) const { return (it != rhs.it); }
			
			iterator& operator+=(difference_type n) { it += n; return *this; }
			iterator& operator-=(difference_type n) { it -= n; return *this; }
			
			pointer operator[](difference_type n) { return (pointer)it[n]; }
			const pointer operator[](difference_type n) const { return (const pointer)it[n]; }
			
			friend iterator operator+(const iterator& i, difference_type n) { return iterator(i.it + n); }
			friend iterator operator+(difference_type n, const iterator& i) { return iterator(i.it + n); }
			friend iterator operator-(const iterator& i, difference_type n) { return iterator(i.it - n); }
			friend iterator operator-(difference_type n, const iterator& i) { return iterator(i.it - n); }
			
			friend difference_type operator-(const iterator& b, const iterator& a) { return (b.it - a.it); }
			
			friend bool operator<(const iterator&a, const iterator& b) { return (a.it < b.it); }
			friend bool operator>(const iterator&a, const iterator& b) { return (a.it > b.it); }
			friend bool operator<=(const iterator&a, const iterator& b) { return (a.it <= b.it); }
			friend bool operator>=(const iterator&a, const iterator& b) { return (a.it >= b.it); }
		};
		
	public:
		/// default constructor
		MixedVector()
		{
		}
		
		/// copy constructor
		MixedVector(const MixedVector<BaseT*>& other)
		: items(other.items), tmap(other.tmap)
		{
		}
		
		/// destroy the vector and all contained items
		~MixedVector()
		{
			clear();
		}
		
		/// swap two MixedVectors
		friend void swap(MixedVector<BaseT*>& a, MixedVector<BaseT*>& b)
		{
			using std::swap;
			
			swap(a.items, b.items);
			swap(a.tmap, b.tmap);
		}
		
		/// assign one MixedVector to another
		MixedVector<BaseT*>& operator=(MixedVector<BaseT*> rhs)
		{
			// copy-and-swap idiom
			swap(*this, rhs);
			return *this;
		}
		
		/// add an object of type SubclassT to the vector
		void push_back(BaseT* obj)
		{
			// make sure this is a subclass of BaseT
			check_subclass(obj);
			
			// add the copy of the object to the list of items
			items.push_back(obj);
			
			// add the copy to type map
			tmap[typeid(*obj)].push_back(obj);
		}
		
		/// pop element from back of vector
		void pop_back()
		{
			// get very last element
			BaseT* back = items.back();
			
			// remove from the vector in the type map
			tmap[typeid(*back)].pop_back();
			
			// remove from items
			items.pop_back();
			
			// delete from memory
			delete back;
		}
		
		/// erase the last element of a specific type
		template<typename SubclassT>
		void pop_back()
		{
			// search for vector belonging to SubclassT
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return;
			
			// pop from vector for type
			found->second.pop_back();
			
			// go through list of all items from back, until item of correct type is found
			for(typename std::vector<BaseT*>::reverse_iterator it = items.rbegin(); it != items.rend(); ++it)
			{
				if(typeid(**it) == typeid(SubclassT))
				{
					// remove the item
					delete(*it);
					items.erase(it);
					return;
				}
			}
		}
		
		/// clear all elements
		void clear()
		{
			while(!empty()) pop_back();
		}
		
		/// clear all elements of a given type
		template<typename SubclassT>
		void clear()
		{
			// remove vector for subclass
			tmap.erase(typeid(SubclassT));
			
			// erase items matching type
			items.erase(std::remove_if(items.begin(), items.end(), is_type<SubclassT>), items.end());
		}
		
		/// direct access to underlying array
		BaseT** data()
		{
			return items.data();
		}
		
		/// direct access to underlying array (const)
		const BaseT** data() const
		{
			return items.data();
		}
		
		/// Checks if element i is of the derived type SubclassT
		template<typename SubclassT>
		bool type(std::size_t i)
		{
			return is_type<SubclassT>(items[i]);
		}
		
		/// Indexing operator for all elements in form of a reference to the base class.
		BaseT* operator[](std::size_t i) const
		{
			return items[i];
		}
		
		/// indexed access
		BaseT* get(std::size_t i)
		{
			return items[i];
		}
		
		/// indexed access (const)
		const BaseT* get(std::size_t i) const
		{
			return items[i];
		}
		
		/// indexed access for type SubclassT
		template<typename SubclassT>
		SubclassT* get(std::size_t i)
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (SubclassT*)found->second[i];
		}
		
		/// indexed access for type SubclassT (const)
		template<typename SubclassT>
		const SubclassT* get(std::size_t i) const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (const SubclassT*)found->second[i];
		}
		
		/// indexed access with bounds checking
		BaseT* at(std::size_t i)
		{
			return items.at(i);
		}
		
		/// indexed access with bounds checking (const)
		const BaseT* at(std::size_t i) const
		{
			return items.at(i);
		}
		
		/** \brief Templated indexing operator for elements of a specific derived class.
		 * The index runs 0 ... size<SubclassT>() - 1
		 * Does bounds checking.
		 */
		template<typename SubclassT>
		SubclassT* at(std::size_t i)
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (SubclassT*)found->second.at(i);
		}
		
		/// Templated indexing operator for elements of a specific derived class (const).
		template<typename SubclassT>
		const SubclassT* at(std::size_t i) const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return (const SubclassT*)found->second.at(i);
		}
		
		/// number of all elements
		std::size_t size() const
		{
			return items.size();
		}
		
		/// number of elements of a specific type
		template<typename SubclassT>
		std::size_t size() const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return 0;
			return found->second.size();
		}
		
		/// check if vector of items is empty
		bool empty() const
		{
			return items.empty();
		}
		
		/// check if vector of items of type SubclassT is empty
		template<typename SubclassT>
		bool empty() const
		{
			typename tmap_t::const_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return true;
			
			return found->second.empty();
		}
		
		/// get iterator to first of all items
		iterator<> begin()
		{
			return iterator<>(items.begin());
		}
		
		/// get iterator to last of all items
		iterator<> end()
		{
			return iterator<>(items.end());
		}
		
		/// get iterator to first of items of type SubclassT
		template<typename SubclassT>
		iterator<SubclassT> begin()
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return iterator<SubclassT>(items.end());
			return iterator<SubclassT>(found->second.begin());
		}
		
		/// get iterator to last of items of type SubclassT
		template<typename SubclassT>
		iterator<SubclassT> end()
		{
			typename tmap_t::iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return iterator<SubclassT>(items.end());
			return iterator<SubclassT>(found->second.end());
		}
		
	private:
		inline void check_subclass(const BaseT*) {}
		
		template<typename SubclassT>
		inline static bool is_type(const BaseT* obj) { return (typeid(*obj) == typeid(SubclassT)); }
		
		typedef std::vector<BaseT*> base_container;
		
		base_container items;
		
		typedef std::map<detail::type_index, base_container> tmap_t;
		tmap_t tmap;
	};
	
	template<class BaseT>
	std::size_t lower_bound(std::vector<BaseT*>& items, double target){
		std::size_t ju,jm,jl;
		
		jl=0;
		ju=items.size()-1;
		while (ju-jl > 1) {
			jm=(ju+jl) >> 1;
			if(items[jm]->compare(target))
				jl=jm;
			else
				ju=jm;
		}
		return jl;
	}
	
	template<typename Container>
	void delete_container(Container& c) { while(!c.empty()) delete c.back(), c.pop_back(); }
  
  /** \brief Class for calculating the Hilbert curve distance in two dimensions
   *
   *  The Hilbert Curve maps two dimensional positions on a grid into distance in one dimension.
   *  The distance d is generally close for points that are near eachother.  This allows one to
   *  uniquely order points in a two dimensional space.  It is useful for matching objects quickly
   *  if they are sorted by thier Hilbert distance.
   *
   */
  class HilbertCurve{
  public:
    HilbertCurve(
                 double x_min       /// x coordinate of lower left corner
                 ,double y_min      /// y coordinate of lower left corner
                 ,double my_range   /// range in which points will be distributed
                 ,double smallsize  /// smallest distance of separation that must give unique Hilbert distances
                 ):
    range(my_range)
    {
      xo[0] = x_min;
      xo[1] = y_min;
      n = (long)(range/smallsize+1);
    }
    
    long xy2d (long x, long y);
    long xy2d (double x, double y);
    void d2xy(long d, long *x, long *y);
    void d2xy(long d, double *x, double *y);
    
  private:
    long n;
    double xo[2],range;
    void rot(long s,long *x, long *y, long rx, long ry);
  };

	template<typename T>
	T between(const T& x, const T& l, const T& u)
	{
		return std::max(l, std::min(u, x));
	}
	
  /// This is a class for generating random numbers. It is actually just a rapper for some std random classes.
  class RandomNumbers{
  public:
    
    RandomNumbers(unsigned int seed);
    ~RandomNumbers(void);
    
    double operator()(void);
    /// Normally (Gaussian) distributed random numbers with mean 0 and standard deviation 1
    double gauss(){return norm_dist(rand_gen);}
  private:
     
    std::normal_distribution<> norm_dist;
    std::mt19937 rand_gen;
  };

}

#endif
