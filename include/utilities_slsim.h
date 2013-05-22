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

namespace Utilities
{
	// this is not for the user
	namespace detail
	{
		class type_index
		{
		public:
			type_index(const std::type_info& type) : t(type) {}
			inline bool operator<(const type_index& rhs) const { return t.before(rhs.t); }
			
		private:
			const std::type_info& t;
		};
	}
	
	//typedef double PosType;
	
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
	 * cout << "Number of Uniform Sources " << mvector.size<SourceUniform>() << "   Number of Gausssian Sources "
	 *		<< mvector.size<SourceGaussian>() << endl;
     *
     * // change derived class attribute
	 * mvector.at<SourceGaussian>(0).source_gauss_r2 = 0.5;
	 * cout << "A base class attribute " << mvector[2].getTotalFlux()
	 *		<< " A derived class attribute  "
	 *		<< mvector.at<SourceGaussian>(0).source_gauss_r2
	 *		<< "  " << mvector.at<SourceGaussian>(1).source_gauss_r2
	 *		<< endl;
	 *
	 *</pre>
	 */
	template<class BaseT>
	class MixedVector
	{
	public:
		/// default constructor
		MixedVector()
		{
		}
		
		/// copy the vector and all contained elements
		MixedVector(const MixedVector<BaseT>& other)
		{
			// needs to clone all the contained base_holder items
			for(typename std::vector<base_holder*>::const_iterator it = other.items.begin(); it != other.items.end(); ++it)
			{
				// clone the object
				base_holder* hold = (*it)->clone();
				
				// add to own items
				items.push_back(hold);
				
				// add to own typemap
				tmap[hold->type()].push_back(hold);
			}
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
		void push_back_ref(const SubclassT& obj)
		{
			// make sure this is a subclass of BaseT
			check_subclass(obj);
			
			// copy the object
			base_holder* hold = new holder<SubclassT>(obj);
			
			// add the copy of the object to the list of items
			items.push_back(hold);
			
			// add the copy to type map
			tmap[typeid(SubclassT)].push_back(hold);
		}
		
		/// pop element from back of vector
		void pop_back()
		{
			// get very last element
			base_holder* back = items.back();
			
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return;
			
			// pop from vector for type
			found->second.pop_back();
			
			// go through list of all items from back, until item of correct type is found
			for(typename std::vector<base_holder*>::reverse_iterator it = items.rbegin(); it != items.rend(); ++it)
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
		
		/// Indexing operator for all elements in form of a reference to the base class.
		BaseT& operator[](std::size_t i) const
		{
			return *(items[i]->to_ptr());
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
			// remove vector for subclass
			tmap.erase(typeid(SubclassT));
			
			// erase items matching type
			items.erase(std::remove_if(items.begin(), items.end(), is_type<SubclassT>), items.end());
		}
		
		/// Checks if element i is of the derived type SubclassT
		template<typename SubclassT>
		bool CheckType(std::size_t i)
		{
			return is_type<SubclassT>(items[i]);
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return *found->second.at(i);
		}
		
		/// Templated indexing operator for elements of a specific derived class (const).
		template<typename SubclassT>
		const SubclassT& at(std::size_t i) const
		{
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return *found->second.at(i);
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return true;
			
			return found->second.empty();
		}
		
		/// get vector of all items
		std::vector<BaseT*> vector() const
		{
			std::vector<BaseT*> results;
			std::transform(items.begin(), items.end(), std::back_inserter(results), to_type<BaseT>);
			return results;
		}
		
		/// get vector of all items of type SubclassT
		template<typename SubclassT>
		std::vector<SubclassT*> vector() const
		{
			std::vector<SubclassT*> results;
			
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return results;
			
			results.reserve(found->second.size());
			std::transform(found->second.begin(), found->second.end(), std::back_inserter(results), to_type<SubclassT>);
			
			return results;
		}
		
	private:
		class base_holder
		{
		public:
			virtual ~base_holder() {}
			
			virtual base_holder* clone() const = 0;
			virtual const std::type_info& type() const = 0;
			virtual BaseT* to_ptr() = 0;
			virtual const BaseT* to_ptr() const = 0;
		};
		
		template<typename T>
		class holder : public base_holder
		{
		public:
			holder(const T& held) { ptr = new T(held); }
			~holder() { delete ptr; }
			
			holder* clone() const { return new holder<T>(*ptr); }
			const std::type_info& type() const { return typeid(T); }
			T* to_ptr() { return ptr; }
			const T* to_ptr() const { return ptr; }
			
		private:
			T* ptr;
		};
		
		inline static void check_subclass(const BaseT&) {}
		
		template<typename SubclassT>
		inline static bool is_type(const base_holder* obj) { return (obj->type() == typeid(SubclassT)); }
		
		template<typename SubclassT>
		inline static SubclassT* to_type(base_holder* obj) { return (SubclassT*)obj->to_ptr(); }
		
		std::vector<base_holder*> items;
		std::map<detail::type_index, std::vector<base_holder*> > tmap;
		
		typedef typename std::map<detail::type_index, std::vector<base_holder*> >::const_iterator tmap_iterator;
	};
	
	/// A MixedVector for pointers.
	template<class BaseT>
	class MixedVector<BaseT*>
	{
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
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
		
		/// Indexing operator for all elements in form of a reference to the base class.
		BaseT* operator[](std::size_t i) const
		{
			return items[i];
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
		
		/// Checks if element i is of the derived type SubclassT
		template<typename SubclassT>
		bool CheckType(std::size_t i)
		{
			return is_type<SubclassT>(items[i]);
		}
		
		/// pointer to first element of items
		BaseT** data()
		{
			return items.data();
		}

		/// pointer to first element of items
		BaseT** ptr(unsigned long j)
		{
			return &items[j];
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return found->second.at(i);
		}
		
		/// Templated indexing operator for elements of a specific derived class (const).
		template<typename SubclassT>
		const SubclassT* at(std::size_t i) const
		{
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				throw std::out_of_range(std::string() + "type " + typeid(SubclassT).name() + " not in vector");
			return found->second.at(i);
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
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
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return true;
			
			return found->second.empty();
		}
		
		/// get vector of all items
		std::vector<BaseT*> vector() const
		{
			return items;
		}
		
		unsigned long lower_bound(double target){
			unsigned long ju,jm,jl;

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

		/// get vector of all items of type SubclassT
		template<typename SubclassT>
		std::vector<SubclassT*> vector() const
		{
			std::vector<SubclassT*> results;
			
			tmap_iterator found = tmap.find(typeid(SubclassT));
			if(found == tmap.end())
				return results;
			
			results.reserve(found->second.size());
			std::transform(found->second.begin(), found->second.end(), std::back_inserter(results), to_type<SubclassT>);
			
			return results;
		}
		
	private:
		inline void check_subclass(const BaseT*) {}
		
		template<typename SubclassT>
		inline static bool is_type(const BaseT* obj) { return (typeid(*obj) == typeid(SubclassT)); }
		
		template<typename SubclassT>
		inline static SubclassT* to_type(BaseT* obj) { return (SubclassT*)obj; }
		
		std::vector<BaseT*> items;
		std::map<detail::type_index, std::vector<BaseT*> > tmap;
		
		typedef typename std::map<detail::type_index, std::vector<BaseT*> >::const_iterator tmap_iterator;
		typedef typename std::vector<BaseT*>::iterator item_iterator;
	};


	template<class BaseT>
	unsigned long lower_bound(std::vector<BaseT*>& items, double target){
		unsigned long ju,jm,jl;

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
}
#endif
