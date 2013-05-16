#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "standard.h"
#include <typeinfo>

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
	void free_Matrix(T **matrix, long rows, long cols){
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
	 * Utilities::MixedVector<Source,SourceGaussian,SourceUniform> mvector;
     *
	 * mvector.push(source1);
	 * mvector.push(source1);
	 * mvector.push(source2);
	 * mvector.push(source2);
	 * cout << "Number of Uniform Sources " << mvector.size<SourceUniform>() << "   Number of Gausssian Sources "
	 *		<< mvector.size<SourceGaussian>() << endl;
     *
     * // change derived class attribute
	 * mvector.index<SourceGaussian>(0).source_gauss_r2 = 0.5;
	 * cout << "A base class attribute " << mvector[2].getTotalFlux()
	 *		<< " A derived class attribute  "
	 *		<< mvector.index<SourceGaussian>(0).source_gauss_r2
	 *		<< "  " << mvector.index<SourceGaussian>(1).source_gauss_r2
	 *		<< endl;
	 *
	 *<\pre>
	 */
	template <class BaseClass,class DClass1,class DClass2>
	class MixedVector{
	public:
		/// add an object of type DClass1 to the vector
		void push(const DClass1 &obj){
			v_class1.push_back(obj);
		}
		/// add an object of type DClass2 to the vector
		void push(const DClass2 &obj){
			v_class2.push_back(obj);
		}
		/// erase element, all DClass2 elements are erased before DClass1 objects
		void pop_back(){
			if(v_class2.size() > 0) v_class2.pop_back();
			else v_class1.pop_back();
		}

		/// Indexing operator for all elements in form of a reference to the base class.
		BaseClass & operator[](size_t in){
			if(in < v_class1.size()) return v_class1[in];
			return v_class2[in - v_class1.size()];
		}

		/// clear all elements
		void clear(){
			v_class1.clear();
			v_class2.clear();
		}

		/// erase the last element of a specific type
		template <typename T>
		void pop(){
			if(typeid(T) == typeid(DClass1)) v_class1.pop_back();
			if(typeid(T) == typeid(DClass2)) v_class2.pop_back();
		}


		template <typename T>
		T & operator[](size_t in){
			if(typeid(T) == typeid(DClass1)) return dynamic_cast<T &> (v_class1[in]);
			if(typeid(T) == typeid(DClass2)) return dynamic_cast<T &> (v_class2[in]);

			ERROR_MESSAGE();
			std::cout << "MixedVector::[] accessed with invalid type." << std::endl;
			exit(1);
		}


		/** \brief Templated indexing operator for elements of a specific derived class.
		 * The index runs 0 ... size<T>() - 1
		 */
		template <typename T>
		T & index(size_t i){

			if(typeid(T) == typeid(DClass1)) return dynamic_cast<T &>( v_class1[i]);
			if(typeid(T) == typeid(DClass2)) return dynamic_cast<T &>( v_class2[i]);

			ERROR_MESSAGE();
			std::cout << "MixedVector::index<type>() accessed with invalid type." << std::endl;
			exit(1);
		}

		/// get a pointer to a specific derived class element, index runs 0 ... size<DClass1>()-1
		void get(DClass1 *obj,size_t index){
			obj = &(v_class1[index]);
		}
		/// get a pointer to a specific derived class element, index runs 0 ... size<DClass2>()-1
		void get(DClass2 *obj,size_t index){
			obj = &(v_class2[index]);
		}

		/// number of elements of a specific type
		template <typename T>
		size_t size(){
			if(typeid(T) == typeid(DClass1)) return v_class1.size();
			if(typeid(T) == typeid(DClass2)) return v_class2.size();
			return 0;
		}

		/// number of all elements
		size_t size(){
			return v_class1.size() + v_class2.size();
		}

	private:
		std::vector<DClass1> v_class1;
		std::vector<DClass2> v_class2;
	};

}
#endif
