#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "standard.h"
#include <typeinfo>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <iterator>

namespace Utilities{
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
		/// destroy the vector and all contained items
		~MixedVector()
		{
			clear();
		}
		
		/// add an object of type SubclassT to the vector
		template<typename SubclassT>
		void push_back_ref(const SubclassT& obj)
		{
			// make sure this is a subclass of BaseT
			check_type(obj);
			
			// copy the object
			BaseT* copy = new SubclassT(obj);
			
			// add the copy of the object to the list of items
			items.push_back(copy);
			
			// add the copy to type map
			tmap[typeid(SubclassT)].push_back(copy);
		}
		
		/// add an object of type SubclassT to the vector
		template<typename SubclassT>
		void push_back(SubclassT* obj)
		{
			// make sure this is a subclass of BaseT
			//check_type(obj);

			// add the copy of the object to the list of items
			items.push_back(obj);

			// add the copy to type map
			tmap[typeid(SubclassT)].push_back(obj);
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
		BaseT& operator[](std::size_t i) const
		{
			return *items[i];
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
			return items;
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
		class type_index
		{
		public:
			type_index(const std::type_info& type) : t(type) {}
			inline bool operator<(const type_index& rhs) const { return t.before(rhs.t); }
			
		private:
			const std::type_info& t;
		};
		
		template<typename SubclassT>
		inline static bool is_type(const BaseT* obj) { return (typeid(*obj) == typeid(SubclassT)); }
		
		template<typename SubclassT>
		inline static SubclassT* to_type(BaseT* obj) { return (SubclassT*)obj; }
		
		inline static void check_type(const BaseT&) {}
		
		std::vector<BaseT*> items;
		std::map<type_index, std::vector<BaseT*> > tmap;
		
		typedef typename std::map<type_index, std::vector<BaseT*> >::const_iterator tmap_iterator;
	};
	
}
#endif
