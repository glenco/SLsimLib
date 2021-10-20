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
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#endif

namespace Utilities
{
/// convert a string to a numerical value of various types
template<class T>
T to_numeric(const std::string &str) {
  return std::stoi(str);
};
template<>
inline long to_numeric<long>(const std::string &str) {
  return std::stol(str);
};
template<>
inline int to_numeric<int>(const std::string &str) {
  return std::stoi(str);
};
template<>
inline float to_numeric<float>(const std::string &str) {
  return std::stof(str);
};
template<>
inline double to_numeric<double>(const std::string &str) {
  return std::stod(str);
};
//********************************************************

template <typename T>
bool AlwaysTrue(T t){return true;}

template <typename T>
bool AlwaysFalse(T t){return false;}

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

PosType **PosTypeMatrix(size_t rows, size_t cols);
void free_PosTypeMatrix(PosType **matrix, size_t rows, size_t cols);
PosType **PosTypeMatrix(long rows1,long rows2, long cols1, long cols2);
void free_PosTypeMatrix(PosType **matrix, long rows1,long rows2, long cols1, long cols2);

/// 3 dimensional matrix, fixed size
template<typename T>
class D3Matrix{
public:
  D3Matrix(size_t xsize,size_t ysize,size_t zsize)
  :xn(xsize),yn(ysize),zn(zsize){
    array = new T[xn*yn*zn];
  }
  ~D3Matrix(){
    delete[] array;
  }
  
  T& operator()(size_t i,size_t j,size_t k){
    return array[i + j*xn + k*xn*yn];
  }
  
  T& operator()(size_t m){
    return array[m];
  }
  
  size_t xindex(size_t m){
    return m % (xn);
  }
  size_t yindex(size_t m){
    return (m % (xn*yn) ) / xn;
  }
  
  size_t zindex(size_t m){
    return m / (xn*yn);
  }
  
private:
  size_t xn;
  size_t yn;
  size_t zn;
  T *array;
};


/// 2 dimensional matrix
template<typename T>
class D2Matrix{
  
public:
  D2Matrix(size_t xsize,size_t ysize)
  :xn(xsize),yn(ysize){
    array = new T[xn*yn];
  }
  ~D2Matrix(){
    delete[] array;
  }
  
  T& operator()(size_t i,size_t j){
    return array[i + j*xn];
  }
  
  T& operator()(size_t m){
    return array[m];
  }
  
  size_t xindex(size_t m){
    return m % xn;
  }
  size_t yindex(size_t m){
    return m / xn;
  }
  
  
private:
  size_t xn;
  size_t yn;
  T *array;
};

/** \brief Symetric matrix
 
 This is a class to represent symmetric maticies so that they can be accessed as
 normal but takes up n(n+1)/2 in memory.
 */
template <typename T>
class SymmetricMatrix{
  std::vector<T> v;
  int n;
  int m;
public:
  SymmetricMatrix(size_t n):n(n){
    v.resize(n*(n+1)/2);
    m = 2*n-1;
  }
  T& operator()(int i,int j){
    //long k = j + (2*n-1-i)*i/2 ;
    //size_t k = (i <= j ) ? j + (2*n-1-i)*i/2 : i + (2*n-1-j)*j/2;
    size_t k = (i <= j ) ? j + (m-i)*i/2 : i + (m-j)*j/2;
    //assert(k>=0);
    //assert(k < v.size());
    return v[ k ];
  }
  T& operator[](size_t k){
    return v[ k ];
  }
  
  int size(){return n;}
  
  /// convertion from 2d to 1d index
  size_t oned_index(int i,int j){
    return (i <= j ) ? j + (m-i)*i/2 : i + (m-j)*j/2;
  }
};


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
    typedef const ValueT* const_pointer;
    typedef ValueT& reference;
    typedef const ValueT& const_reference;
    typedef typename base_iterator::difference_type difference_type;
    typedef typename base_iterator::iterator_category iterator_category;
    
    iterator() {}
    iterator(base_iterator i) : it(i) {}
    iterator(const iterator& other) : it(other.it) {}
    
    iterator& operator=(const iterator& rhs) { it = rhs.it; return *this; }
    
    reference operator*() { return (reference)(**it); }
    
    //const reference operator*() const { return (const reference)(**it); }
    
    pointer operator->() { return (pointer)(*it); }
    const_pointer operator->() const { return (const_pointer)(*it); }
    
    iterator& operator++() { ++it; return *this; }
    iterator operator++(int) { iterator tmp(*this); ++it; return tmp; }
    iterator& operator--() { --it; return *this; }
    iterator operator--(int) { iterator tmp(*this); --it; return tmp; }
    
    bool operator==(const iterator& rhs) const { return (it == rhs.it); }
    bool operator!=(const iterator& rhs) const { return (it != rhs.it); }
    
    iterator& operator+=(difference_type n) { it += n; return *this; }
    iterator& operator-=(difference_type n) { it -= n; return *this; }
    
    reference operator[](difference_type n) { return (reference)*it[n]; }
    
    //const reference operator[](difference_type n) const { return (const reference)*it[n]; }
    
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
  
  /// back of list
  BaseT& back()
  {
    return items.back();
  }
  
  /// back of list
  const BaseT& back() const
  {
    return items.back();
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
std::size_t lower_bound(std::vector<BaseT*>& items, PosType target){
  std::size_t ju,jm,jl;
  
  if(items.size() == 0) return 0;
  
  jl=0;
  ju=items.size()-1;
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if(items[jm]->compareZ(target))
      jl=jm;
    else
      ju=jm;
  }
  return jl;
}

/// delete the objects that are pointed to in a container of pointers
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
               PosType x_min       /// x coordinate of lower left corner
               ,PosType y_min      /// y coordinate of lower left corner
               ,PosType my_range   /// range in which points will be distributed
               ,PosType smallsize  /// smallest distance of separation that must give unique Hilbert distances
  ):
  range(my_range)
  {
    xo[0] = x_min;
    xo[1] = y_min;
    n = (int)(range/smallsize+1);
  }
  
  int xy2d (int x, int y);
  int xy2d (PosType x, PosType y);
  void d2xy(int d, int *x, int *y);
  void d2xy(int d, PosType *x, PosType *y);
  
private:
  int n;
  PosType xo[2],range;
  void rot(int s,int *x, int *y, int rx, int ry);
};

template<typename T>
T between(const T& x, const T& l, const T& u)
{
  return std::max(l, std::min(u, x));
}

#ifdef ENABLE_CLANG
/// This is a class for generating random numbers. It is actually just a rapper for some std random classes.
class RandomNumbers{
public:
  
  RandomNumbers(unsigned int seed);
  ~RandomNumbers(void);
  
  PosType operator()(void);
  /// Normally (Gaussian) distributed random numbers with mean 0 and standard deviation 1
  PosType gauss(){return norm_dist(rand_gen);}
private:
  
  std::normal_distribution<> norm_dist;
  std::mt19937 rand_gen;
};
#endif

/**
 * \brief This is a class for generating random numbers. It simplifies and fool proofs initialization and allows for multiple
 *  independent series of numbers.
 *
 * This version is based on NR ran2() and is provided only for backwards reproducibility.  Use RandomNumbers class when possible.
 */
class RandomNumbers_NR{
public:
  
  RandomNumbers_NR(long seed);
  
  PosType operator()(void);
  /// generates a Gaussian distributed number with unit variance by polar Box-Muller transform
  PosType gauss(){
    ++calls;
    if(count){
      do{
        u = 2*ran2() - 1;
        v = 2*ran2() - 1;
        s = u*u +v*v;
      }while( s > 1.0 || s == 0.0);
      
      s = sqrt(-2*log(s)/s);
      count = false;
      return s*u;
    }else{
      count = true;
      return s*v;
    }
  };
  
  int poisson(double lam){
    double L = exp(-lam),p=1;
    int k = 0;
    do{
      ++k;
      p *= operator()();
    }while(p>L);
    
    return k-1;
  }
  
  size_t calls = 0;  /// total number of calls
  
  long getseed(){return firstseed;}
private:
  long idum;
  PosType ran2(void);
  
  int IM1;
  int IM2;
  PosType AM;
  //int IMM1 = (IM1-1);
  int IA1;
  int IA2;
  int IQ1;
  int IQ2;
  int IR1;
  int IR2;
  int NDIV;
  PosType EPS;
  PosType RNMX;
  
  long idum2;
  long iy;
  long iv[32];
  bool count;
  PosType u,v,s;
  long firstseed;
  
};

/// Shuffles a vector into a random order
template <typename T, typename R>
void shuffle(
             std::vector<T> &vec   /// The vector to be shuffled
             ,R &ran               /// a random number generator so that ran() gives a number between 0 and 1
){
  T tmp;
  size_t ran_t;
  if(vec.size() < 2) return;
  for (size_t i = vec.size()-1; i>0; --i) {
    ran_t = (size_t)(ran()*(i+1));
    //swap(tmp,vec[ran_t]);
    std::swap(vec[ran_t],vec[i]);
    //swap(vec[i],tmp);
  }
}

/// Find the indexes that sort a vector in asending order
template <typename T>
void sort_indexes(const std::vector<T> &v     /// the original data that is not changed
                  ,std::vector<size_t> &index /// vector of indexes that if put into v will sort it
) {
  
  // initialise original index locations
  index.resize(v.size());
  for (size_t i = 0; i != index.size(); ++i) index[i] = i;
  
  // sort indexes based on comparing values in v
  std::sort(index.begin(), index.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
}

template <typename T>
void sort_indexes(const T *v     /// the original data that is not changed
                  ,std::vector<size_t> &index /// vector of indexes that if put into v will sort it
                  ,size_t N) {
  
  // initialise original index locations
  index.resize(N);
  
  for (size_t i = 0; i != index.size(); ++i) index[i] = i;
  
  // sort indexes based on comparing values in v
  std::sort(index.begin(), index.end(),
            [v](size_t i1, size_t i2) {return v[i1] < v[i2];});
}

/// Find the indexes that sort a vector in descending order
template <typename T>
void sort_indexes_decending(const std::vector<T> &v     /// the original data that is not changed
                            ,std::vector<size_t> &index /// vector of indexes that if put into v will sort it
) {
  
  // initialise original index locations
  index.resize(v.size());
  for (size_t i = 0; i != index.size(); ++i) index[i] = i;
  
  // sort indexes based on comparing values in v
  std::sort(index.begin(), index.end(),
            [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
}

/** \brief Gives a randomized sequence of numbers from 0 to N-1.
 
 If the sequence is exhausted then it is reshuffled and numbers will
 repeat in randomized order.
 **/
template <typename R>
class ShuffledIndex{
public:
  ShuffledIndex(size_t N,R &ran){
    if(N == 0) throw std::invalid_argument("N=0");
    
    for(size_t i=0 ; i<N ; ++i) index[i] = i;
    index.resize(N);
    Utilities::shuffle(index,ran);
    index_internal = 0;
  }
  
  /// return the index number
  size_t operator*(){return index[index_internal];}
  /// goto the next number
  size_t operator++(int){
    size_t tmp = index[index_internal];
    if(index_internal == index.size()){
      index_internal = 0;
    }else{
      ++index_internal;
    }
    return tmp;
  }
  /// goto the next number
  size_t operator++(){
    if(index_internal == index.size()){
      index_internal = 0;
    }else{
      ++index_internal;
    }
    return index[index_internal];
  }
  
  /// get a new random order
  void reshuffle(R &ran){
    Utilities::shuffle(index,ran);
    index_internal = 0;
  }
  
private:
  std::vector<size_t> index;
  size_t index_internal;
};


// reorders vec according to index p
template <typename T>
void apply_permutation(
                       T *vec,
                       const std::vector<std::size_t>& p)
{
  std::vector<T> copy(p.size());
  
  for(std::size_t i = 0; i < p.size(); ++i){
    copy[i] = vec[i];
  }
  
  for(std::size_t i = 0; i < p.size(); ++i){
    vec[i] = copy[p[i]];
  }
}

template <typename T>
void apply_permutation(
                       std::vector<T>& vec,
                       const std::vector<std::size_t>& p)
{
  apply_permutation(vec.data(),p);
}

#ifdef ENABLE_FFTW
/** \brief Calculates power spectrum from a 2d map or the cross-power spectrum between two 2d maps.
 *
 *   Adaptied from Carlo Giocoli's pl() routine.
 */
void powerspectrum2d(
                     std::valarray<double> const &aa      /// first realspace map to be
                     ,std::valarray<double> const &bb     /// second realspace map, same as aa to get power spectrum
                     ,int nx                       /// number of pixels in x direction
                     ,int ny                       /// number of pixels in y direction
                     ,double boxlx                 /// range of image in x direction
                     ,double boxly                 /// range of image in y direction
                     ,std::vector<double> &ll      /// output multiplot number of bins
                     ,std::vector<double> &Pl      /// output binned power spectrum
                     ,double zeropaddingfactor
                     );
void powerspectrum2d(
                     std::valarray<double> &aa      /// first realspace map to be
                     ,int nx                       /// number of pixels in x direction
                     ,int ny                       /// number of pixels in y direction
                     ,double boxlx                 /// range of image in x direction
                     ,double boxly                 /// range of image in y direction
                     ,std::vector<double> &ll      /// output multiplot number of bins
                     ,std::vector<double> &Pl      /// output binned power spectrum
);
void powerspectrum2dprebin(
                     std::valarray<double> &aa      /// first realspace map to be
                     ,int nx                       /// number of pixels in x direction
                     ,int ny                       /// number of pixels in y direction
                     ,double boxlx                 /// range of image in x direction
                     ,double boxly                 /// range of image in y direction
                     ,const std::vector<double> &ll      /// output multiplot number of bins
                     ,std::vector<double> &Pl      /// output binned power spectrum
                     ,std::vector<double> &llave     /// average value of Fourier node in bins
);

#endif


/** \brief Smooth a 2 dimensional map stored in a valarray with a density dependent kernel.
 
 The smoothing is done by finding the circle around each point whose total pixel values are larger than value.  In the case of a density map made from particles if value = (mass of particle)*(number of neighbours) an approximate N nearest neighbour smoothing is done.
 The
 **/
std::valarray<double> AdaptiveSmooth(const std::valarray<double> &map_in,size_t Nx,size_t Ny,double value);

/** \brief Smooth a 2 dimensional map stored in a valarray with a density dependent kernel.
 
 The smoothing is done by finding the circle around each point whose total pixel values are larger than value.  In the case of a density map made from particles if value = (mass of particle)*(number of neighbours) an approximate N nearest neighbour smoothing is done.
 **/
std::vector<double> AdaptiveSmooth(const std::vector<double> &map_in,size_t Nx,size_t Ny,double value);

/// returns the compiler variable N_THREADS that is maximum number of threads to be used.
int GetNThreads();

/// namespace for input/output utilities
namespace IO{  ///

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

/** \brief Read in data from an ASCII file with two columns
 */
template <class T1,class T2>
void read2columnfile(
                     std::string filename    /// input file name
                     ,std::vector<T1> &x     /// vector that will contain the first column
                     ,std::vector<T2> &y     /// vector that will contain the second column
                     ,std::string delineator = " "  /// specific string the seporates columns, ex. ",", "|", etc.
                     ,int skiplines = 0
                     ,bool verbose = false
                     
                     ){
  
  x.clear();
  y.clear();
  
  std::ifstream file_in(filename.c_str());
  std::string myline;
  std::string space = " ";
  T1 myt1;
  T2 myt2;
  
  std::string strg;
  std::stringstream buffer;
  
  if(!file_in){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  
  std::cout << "Reading caustic information from " << filename << std::endl;
  size_t i=0;
  while(i < skiplines){
    getline(file_in,myline);
    ++i;
  }
  while(file_in.peek() == '#'){
    file_in.ignore(10000,'\n');
    ++i;
  }
  std::cout << "skipped "<< i << " comment lines in " << filename << std::endl;
  
  size_t pos;
  // read in data
  while(getline(file_in,myline)){
    
    if(myline[0] == '#'){
      std::cout << "skipped line " << i << std::endl;
      continue;
    }
    
    pos= myline.find_first_not_of(space);
    myline.erase(0,pos);
    
    
    pos = myline.find(delineator);
    strg.assign(myline,0,pos);
    buffer << strg;
    buffer >> myt1;
    if(verbose) std::cout << myt1 << " ";
    x.push_back(myt1);
    
    myline.erase(0,pos+1);
    pos= myline.find_first_not_of(space);
    myline.erase(0,pos);
    
    strg.clear();
    buffer.clear();
    
    pos = myline.find(space);
    strg.assign(myline,0,pos);
    buffer << strg;
    buffer >> myt2;
    if(verbose)  std::cout << myt2 << std::endl;
    y.push_back(myt2);
    
    strg.clear();
    buffer.clear();
    myline.clear();
    
  }
  std::cout << "Read " << x.size() << " lines from " << filename << std::endl;
}
/** \brief Read in data from an ASCII file with three columns
 */
template <class T1,class T2,class T3>
void read3columnfile(
                     std::string filename    /// input file name
                     ,std::vector<T1> &x     /// vector that will contain the first column
                     ,std::vector<T2> &y     /// vector that will contain the second column
                     ,std::vector<T3> &z     /// vector that will contain the third column
                     ,std::string delineator = " "  /// specific string the seporates columns, ex. ",", "|", etc.
                     ,bool verbose = false
                     
                     ){
  
  
  assert(0); // Untested!!!!
  x.clear();
  y.clear();
  z.clear();
  
  std::ifstream file_in(filename.c_str());
  std::string myline;
  std::string space = " ";
  T1 myt1;
  T2 myt2;
  T3 myt3;
  
  std::string strg;
  std::stringstream buffer;
  
  if(!file_in){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  
  std::cout << "Reading caustic information from " << filename << std::endl;
  size_t i=0;
  while(file_in.peek() == '#'){
    file_in.ignore(10000,'\n');
    ++i;
  }
  std::cout << "skipped "<< i << " comment lines in " << filename << std::endl;
  
  size_t pos;
  // read in data
  while(getline(file_in,myline)){
    
    if(myline[0] == '#'){
      std::cout << "skipped line " << i << std::endl;
      continue;
    }
    
    pos= myline.find_first_not_of(space);
    myline.erase(0,pos);
    
    
    pos = myline.find(delineator);
    strg.assign(myline,0,pos);
    buffer << strg;
    buffer >> myt1;
    if(verbose) std::cout << myt1 << " ";
    x.push_back(myt1);
    
    myline.erase(0,pos+1);
    pos= myline.find_first_not_of(space);
    myline.erase(0,pos);
    
    strg.clear();
    buffer.clear();
    
    // ******************
    
    
    pos = myline.find(space);
    strg.assign(myline,0,pos);
    buffer << strg;
    buffer >> myt2;
    if(verbose)  std::cout << myt2 << std::endl;
    y.push_back(myt2);
    
    myline.erase(0,pos+1);
    pos= myline.find_first_not_of(space);
    myline.erase(0,pos);
    
    strg.clear();
    buffer.clear();
    
    // ******************
    pos = myline.find(space);
    strg.assign(myline,0,pos);
    buffer << strg;
    buffer >> myt3;
    if(verbose)  std::cout << myt3 << std::endl;
    y.push_back(myt3);
    
    strg.clear();
    buffer.clear();
    myline.clear();
    
  }
  std::cout << "Read " << x.size() << " lines from " << filename << std::endl;
}

int NumberOfEntries(const std::string &string,char deliniator);

/// Count the number of columns in a ASCII data file.
int CountColumns(std::string filename  /// name of file
                 ,char comment_char = '#'  /// comment charactor
                 ,char deliniator = ' '    /// deliniator between columns
);

/** \brief Reads the file names in a directory that contain a specific sub string.
 
 */
void ReadFileNames(
                   std::string dir              /// path to directory containing fits files
                   ,const std::string filespec /// string of charactors in file name that are matched. It can be an empty string.
                   ,std::vector<std::string> & filenames  /// output vector of PixelMaps
                   ,bool verbose);

/// check if the directory does not exist
bool check_directory(std::string dir);


/** \brief This function will read in all the numbers from a multi-column
 ,space seporated ASCII data file.
 
 It will skip the comment lines if they are at the head of the file.  The
 number of columns and rows are returned.  The entry at row r and column c will be stored at data[c + column*r].
 
 This function is not particularly fast for large amounts of data.  If the
 number of roaws is large it would be best to use data.reserve() to set the capacity of data large enough that no rellocation of memory occurs.
 */
template <typename T>
void ReadASCII(std::vector<T> &data
               ,std::string filename
               ,int &columns
               ,int &rows
               ,char comment_char = '#'
               ,int skiplines = 0
               ,size_t MaxNrows = std::numeric_limits<size_t>::max()
               ,bool verbose = true){
  
  std::ifstream file(filename);
  // find number of particles
  if (!file.is_open()){
    std::cerr << "file " << filename << " cann't be opened." << std::endl;
    throw std::runtime_error("no file");
  }
  
  
  data.empty();
  
  std::string line;
  columns = 0;
  rows = 0;
  size_t first_data_line;
  
  // skip over first lines
  int i=0;
  while(i < skiplines){
    std::getline(file,line);
    if(!file) break;  // probably EOF
    ++i;
  }
  
  // read comment lines and first data line
  do{
    first_data_line = file.tellg();
    std::getline(file,line);
    if(!file) break;  // probably EOF
  }while(line[0] == comment_char);
  
  columns =  NumberOfEntries(line,' ');
  
  file.seekg(first_data_line);   // move back to first data line
  
  std::copy(std::istream_iterator<T>(file),
            std::istream_iterator<T>(),
            std::back_inserter(data));
  
  rows = data.size()/columns;
  if(verbose){
    std::cout << "Read " << rows << " rows of " << columns << " columns from file " << filename << std::endl;
  }
}

/** \brief Read numerical data from a csv file with a header
 
 It will skip the comment lines if they are at the head of the file.  The
 number of columns and rows are returned.  The entries will be stored at data[column][row].
 
 Comments must only be before the data.  There must be a line with the
 column names after the comments and before the data.
 
 This function is not particularly fast for large amounts of data.  If the
 number of rows is large it would be best to use data.reserve() to set the capacity of data large enough that no rellocation of memory occurs.
 
 * The accept function can be used to limit the amount of data added.  If there is an object, a, used to
 *    make this selection this can be done like [&a](str::vector<T> &v}{return a.itsok(v[3],v[4]);}
 *         where v corresponds to a row in the data file in order.
 
 */

template <typename T>
int ReadCSVnumerical1(std::string filename                              /// file name to be read
                      ,std::vector<std::vector<T> > &data               /// output data
                      ,std::vector<std::string> &column_names           /// list of column names
                       ,size_t MaxNumber = 100000000                     /// maximum number of entries read
                      ,char comment_char = '#'                          /// comment charactor for header
                      ,char deliniator = ','                            /// deliniator between values
                      ,std::string replace = "\\N"                      /// replace this string with zero
                      ,std::function<bool(std::vector<T> &)> accept = [](std::vector<T> &v){return true;}  /// function that determines if a row should be accepted
){
  
  bool verbose = false;
  
  std::ifstream file(filename);
  // find number of particles
  if (!file.is_open()){
    std::cerr << "file " << filename << " cann't be opened." << std::endl;
    throw std::runtime_error("no file");
  }
  
  int count_preamble=0;
  std::string line;
  // read comment lines and first data line
  do{
    std::getline(file,line);
    ++count_preamble;
    if(!file) break;  // probably EOF
  }while(line[0] == comment_char);
  
  // read the names
  std::stringstream          lineStream(line);
  std::string                cell;
  column_names.empty();
  while(std::getline(lineStream,cell, deliniator))
  {
    column_names.push_back(cell);
  }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty())
  {
    column_names.push_back("");
  }
  
  if(verbose){ // print colum names
    int i = 0;
    for(auto st : column_names){
      std::cout << i++ << " " << st << std::endl;
    }
  }
  
  int columns = NumberOfEntries(line,deliniator);
  
  // count number of data line
  size_t number_of_data_lines = 0;
  while(std::getline(file, line) && number_of_data_lines < MaxNumber ) ++number_of_data_lines;
  
  //std::vector<std::vector<T> > tmp_data(columns,std::vector<T>(number_of_data_lines));
  std::vector<std::vector<T> > tmp_data(columns,std::vector<T>(0));
  swap(data,tmp_data);
  std::vector<T> tmp_row(columns);
  
  /// return to first data line
  //file.seekg(std::ios::beg);
  file.clear();
  file.close();
  file.open(filename);
  for(int i=0; i < count_preamble; ++i){
    file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  //std::getline(file,line);
  
  size_t rows = 0;
  //data.resize(columns);
  //for(auto &v : data ) v.empty();
  while(std::getline(file,line) && rows < MaxNumber){
    // read the names
    std::stringstream          lineStream(line);
    std::string                cell;
    
    int i=0;
    while(std::getline(lineStream,cell,deliniator))
    {
      if(cell==replace) cell='0';
      
      /// clean blank spaces
      cell.erase(remove_if(cell.begin(),cell.end(), isspace), cell.end());
      //data[i].push_back(to_numeric<T>(cell));
      //data[i][rows] = to_numeric<T>(cell);
      tmp_row[i] = to_numeric<T>(cell);
      i = (i+1)%columns;
    }
    if(accept(tmp_row)){
      ++rows;
      for(int i=0 ; i < columns ; ++i) data[i].push_back(tmp_row[i]);
    }
  }
  return 1;
}

/**
 Finds ranges
 */

template <typename T>
size_t ReadCSVrange(std::string filename                              /// file name to be read
                      ,std::vector<std::vector<T> > &ranges               /// output data
                      ,std::vector<std::string> &column_names           /// list of column names
                       ,size_t MaxNumber = 100000000                     /// maximum number of entries read
                      ,char comment_char = '#'                          /// comment charactor for header
                      ,char deliniator = ','                            /// deliniator between values
                      ,std::string replace = "\\N"                      /// replace this string with zero
                      ,std::function<bool(std::vector<T> &)> accept = [](std::vector<T> &v){return true;}  /// function that determines if a row should be accepted
){
  
  bool verbose = false;
  
  std::ifstream file(filename);
  // find number of particles
  if (!file.is_open()){
    std::cerr << "file " << filename << " cann't be opened." << std::endl;
    throw std::runtime_error("no file");
  }
  
  int count_preamble=0;
  std::string line;
  // read comment lines and first data line
  do{
    std::getline(file,line);
    ++count_preamble;
    if(!file) break;  // probably EOF
  }while(line[0] == comment_char);
  
  // read the names
  std::stringstream          lineStream(line);
  std::string                cell;
  column_names.empty();
  while(std::getline(lineStream,cell, ','))
  {
    column_names.push_back(cell);
  }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty())
  {
    column_names.push_back("");
  }
  
  if(verbose){ // print colum names
    int i = 0;
    for(auto st : column_names){
      std::cout << i++ << " " << st << std::endl;
    }
  }
  
  int columns = NumberOfEntries(line,deliniator);
  
  std::vector<std::vector<T> > tmp_ranges(columns,std::vector<T>(2));
  swap(ranges,tmp_ranges);
  std::vector<T> tmp_row(columns);
  
  /// return to first data line
  //file.seekg(std::ios::beg);
  file.clear();
  file.close();
  file.open(filename);
  for(int i=0; i < count_preamble; ++i){
    file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  //std::getline(file,line);
  
  size_t rows = 0;
  while(std::getline(file,line) && rows < MaxNumber){
    // read the names
    std::stringstream          lineStream(line);
    std::string                cell;
    
    int i=0;
    while(std::getline(lineStream,cell,deliniator))
    {
      if(cell==replace) cell='0';
      
      /// clean blank spaces
      cell.erase(remove_if(cell.begin(),cell.end(), isspace), cell.end());
      tmp_row[i] = to_numeric<T>(cell);
      i = (i+1)%columns;
    }
    if(accept(tmp_row)){
      if(rows==0){
        for(int j=0 ; j < columns ; ++j){
          ranges[j][0] = tmp_row[j];
          ranges[j][1] = tmp_row[j];

        }
      }else{
        for(int j=0 ; j < columns ; ++j){
          if(ranges[j][0] > tmp_row[j]) ranges[j][0] = tmp_row[j];
          if(ranges[j][1] < tmp_row[j]) ranges[j][1] = tmp_row[j];
        }
      }
      ++rows;
    }
  }
  return rows;
}

/** \brief Read numerical data from a csv file with a header
 
 Same as ReadCSVnumerical1 except the order of the data storage is
 reversed data[row][column].
 
 It will skip the comment lines if they are at the head of the file.  The
 number of columns and rows are returned.
 
 Comments must only be before the data.  Then if header==true there
 should be a line of column names.  If header!=true there are non column names.
 
 This function is not particularly fast for large amounts of data.  If the
 number of rows is large it would be best to use data.reserve() to set the capacity of data large enough that no rellocation of memory occurs.
 */

template <typename T>
int ReadCSVnumerical2(std::string filename   /// file name to be read
                      ,std::vector<std::vector<T> > &data  /// output data
                      ,std::vector<std::string> &column_names /// list of column names
                      ,size_t Nmax = 1000000
                      ,char comment_char = '#'  /// comment charactor for header
                      ,char deliniator = ','    /// deliniator between values
                      ,bool header = true    /// false if there are no column names
                      ,std::string reject = ""  /// reject lines with this entry after striping black space
){
  
  
  std::ifstream file(filename);
  
  if (!file.is_open()){
    std::cerr << "file " << filename << " cann't be opened." << std::endl;
    throw std::runtime_error("no file");
  }
  std::string line;
  // read comment lines and first data line
  do{
    std::getline(file,line);
    if(!file) break;  // probably EOF
  }while(line[0] == comment_char);
  
  // read the names
  std::stringstream          lineStream(line);
  std::string                cell;
  column_names.empty();
  long ii = 0;
  
  if(header){
    while(std::getline(lineStream,cell,deliniator))
    {
      column_names.push_back(cell);
    }
    
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
      column_names.push_back("");
    }
  }else{
    
    while(std::getline(lineStream,cell,deliniator))
    {
      column_names.push_back(cell);
    }
    
    data.emplace_back(column_names.size());
    int i=0;
    for(auto a : column_names){
      /// clean blank spaces
      cell.erase(remove_if(cell.begin(),cell.end(), isspace), cell.end());
      data.back()[i++] = to_numeric<T>(a);
    }
    ++ii;
  }
  
  int columns = NumberOfEntries(line,deliniator);
  
  for(auto &v : data ) v.empty();
  
  while(std::getline(file,line) && ii < Nmax){
    
    data.emplace_back(columns);
    
    // read the numbers
    std::stringstream          lineStream(line);
    std::string                cell;
    
    int i=0;
    while(std::getline(lineStream,cell, deliniator))
    {
      assert(i < columns);
      cell.erase(remove_if(cell.begin(),cell.end(), isspace), cell.end());
      if(cell == reject){
        data.pop_back();
        --ii;
        break;
      }else{
        data.back()[i] = to_numeric<T>(cell);
        ++i;
      }
    }
    ++ii;
  }
  
  return 1;
}

/** \brief write a CSV data file for some data vectors
 
 example:
 <p>
 std::vector<std::string> header = {"alpha","kappa","gamma"};
 
 std::vector<double> v1 = {1,2,3};
 std::vector<double> v2 = {1.1,2.1,3.1};
 std::vector<double> v3 = {3,4,5};
 
 std::vector<std::vector<double> *> data;
 data.push_back(&v1);
 data.push_back(&v2);
 data.push_back(&v3);
 writeCSV(filename,header,data);
 </p>
 **/

template<typename T>
void writeCSV(const std::string filename              /// output file path/name
              ,const std::vector<std::string> header  /// column labels
              ,std::vector<T *> &data                 /// objects must have operator []
){
  
  std::ofstream s(filename + ".csv");
  
  int ncol = header.size();
  assert(ncol == data.size() );
  for(int i = 0 ; i < header.size()-1 ; ++i){
    s << header[i] << ",";
  }
  s << header.back() << std::endl;
  
  size_t nrow = data[0]->size();
  for(auto v : data) assert(nrow == v->size() );
  
  for(size_t j = 0 ; j < nrow ; ++j){
    for(int i = 0 ; i < ncol-1 ; ++i){
      s << data[i]->operator[](j) << ",";
    }
    s << data.back()->operator[](j) << "\n";
  }
}

/*** \brief A class for reading and then looking up objects from a CSV catalog.
 
 The constructor will read in the whole catalog and sort it into Nxbins X-bins.  Each
 X-bin is then sorted by Y.
 
 The find(x,y) function will set the current value to a galaxy with a x within the
 x-bin of x and a y close to y.  The other information in the row of the csv file for
 this entry can then be read.
 
 */
class XYcsvLookUp{
public:
  XYcsvLookUp(std::string datafile   /// input catalog file in csv format
              ,std::string Xlabel    /// label in catalog for X variable
              ,std::string Ylabel    /// label in catalog for Y variable
              ,int Nxbins            /// number of X bins
              ,size_t MaxNumber = 1000000 /// maximum number of entries read
              ,bool verbose = false);
  XYcsvLookUp(
              std::string datafile   /// input catalog file in csv format
              ,std::string Xlabel    /// label in catalog for X variable
              ,std::string Ylabel    /// label in catalog for Y variable
              ,std::vector<double> Xbins  /// minimum X value in each bin
              ,size_t MaxNumber = 1000000 /// maximum number of entries read
              ,bool verbose = false);
  
  /// find a galaxy with a redshift and log(halo mass) near x and logm, return a vector of its properties
  std::vector<double> find(double x,double y);
  
  /// returns the current galxay's property by label
  double operator[](std::string label) const;
  
  /// returns the ith entry for the current galaxy
  double operator[](int i) const{
    return (*current)[i];
  }
  /// returns the redshift of the current galaxy
  double getX() const {
    return (*current)[Xindex];
  }
  /// returns the log(mass) of the current halo
  double getY() const {
    return (*current)[Yindex];
  }
  /// returns a vector of the entries for the current galaxy
  std::vector<double> record() const{
    return *current;
  }
  /// returns a vector of the entries for the current galaxy
  std::vector<double> operator*() const {
    return *current;
  }
  
  void operator++(){
    if(current != data.end() ) ++current;
  }
  void operator--(){
    if(current != data.begin() ) --current;
  }
  /// returns labels of the columns from the data file
  std::vector<std::string> labels() const{
    return column_names;
  }
  
  double Xmin() const {return xmin;}
  double Xmax() const {return xmax;}
  /// returns minimum Y value in x-bin
  double Ymin(double x) const;
  /// returns maximum Y value in x-bin
  double Ymax(double x) const;
  
  void printcurrent(){
    int i = 0;
    for(auto label : column_names){
      std::cout << label << " : " << (*current)[i++] << std::endl;
    }
  }
  
private:
  int Xindex;
  int Yindex;
  double xmax;
  double xmin;
  std::vector<std::vector<double> >::iterator current;
  std::vector<std::vector<std::vector<double> >::iterator> borders;
  std::vector<std::vector<double> > data;
  std::vector<std::string> column_names;
  std::vector<double> Xborders;
  size_t NinXbins;
  std::string filename;
};
} // Utilities::IO

/// split string into vector of seporate strings that were seporated by
void splitstring(std::string &line,std::vector<std::string> &vec,const std::string &delimiter);

/** \brief class for impoting data from a csv file and allowing label string lookup like a data frame.
 *
 * The accept function can be used to limit the amount of data added.  If there is an object, a, used to
 *    make this selection this can be done like [&a](str::vector<T> &v}{return a.itsok(v[3],v[4]);}
 *         where v corresponds to a row in the data file in order.
 */
template< typename T>
class DataFrame{
public:
  DataFrame(std::string datafile   /// input catalog file in csv format
            ,size_t MaxNumber = 1000000 /// maximum number of entries read
            ,char comment_char = '#'  /// comment charactor for header
            ,char deliniator = ','    /// deliniator between values
            ,std::string replace = "\\N"    /// replace this string with zeros
            ,std::function<bool(std::vector<T> &)> accept = [](std::vector<T> &v){return true;}  /// function that determines if a row should be accepted
  ):filename(datafile){
    input(datafile,MaxNumber,comment_char,deliniator,replace,accept);
  };
 
  DataFrame(){};
  
  void input(std::string datafile   /// input catalog file in csv format
            ,size_t MaxNumber = 1000000 /// maximum number of entries read
            ,char comment_char = '#'  /// comment charactor for header
            ,char deliniator = ','    /// deliniator between values
            ,std::string replace = "\\N"    /// replace this string with zeros
            ,std::function<bool(std::vector<T> &)> accept = [](std::vector<T> &v){return true;}  /// function that determines if a row should be accepted
  ){
    filename = datafile;
    Utilities::IO::ReadCSVnumerical1(datafile,data,column_names,MaxNumber,'#',',',replace,accept);
    
    for(int i=0 ; i<column_names.size() ; ++i) datamap[column_names[i]] = i;
  };

  /// remove a column
  void pop(std::string colname){
    
    int i=datamap[colname];
 
    swap(column_names[i],column_names.back());
    column_names.pop_back();
    swap(data[i],data.back());
    data.pop_back();
    
    datamap.empty();
    for(int i=0 ; i<column_names.size() ; ++i) datamap[column_names[i]] = i;
  }
  
  /// returns column by name
  std::vector<T>& operator[](const std::string &label){
    if(datamap.find(label) == datamap.end()){
      std::cerr << "No label - " << label << " - in " << filename <<std::endl;
      throw std::invalid_argument("no label");
    }
    return data[datamap[label]];
    for(auto c : column_names ) std::cout << c << " ";
    std::cout << std::endl;
    throw std::invalid_argument(label + " was not one of the columns of the galaxy data file :" + filename);
  };
  
  /// add a column to the data frame.  This does a copy.
  void add_column(const std::string &name,const std::vector<T> &vec){
    if(vec.size() != data[0].size()){
      std::cerr << "Added column to DataFrame needs to have the same size." << std::endl;
      throw std::invalid_argument("wrong size");
    }
    column_names.push_back(name);
    data.push_back(vec);
    data[name] = data.size()-1;
  };
  
  /// returns column by number
  std::vector<T>& operator[](int i){
    return data[i];
  };
  
  /// returns labels of the columns from the data file
  std::vector<std::string> labels() const{
    return column_names;
  };
  
  // sort by one of the columns
  void sortby(std::string name){
    std::vector<size_t> index(data[0].size());
    size_t N = index.size();
    for(size_t i=0 ; i<N ; ++i) index[i] = i;
    
    sort_indexes(data[datamap[name]],index);
    std::vector<T> tmp_v(N);
    
    for(size_t j=0 ; j<data.size() ; ++j){
      for(size_t i=0 ; i<N ; ++i){
        tmp_v[i] = data[j][index[i]];
      }
      swap(data[j],tmp_v);
    }
  }
  
  // shuffle order of rows
  void shuffle(Utilities::RandomNumbers_NR &ran){
    std::vector<size_t> index(data[0].size());
    size_t N = index.size();
    for(size_t i=0 ; i<N ; ++i) index[i] = i;
    Utilities::shuffle(index, ran);
    
    std::vector<T> tmp_v(N);
    
    for(size_t j=0 ; j<data.size() ; ++j){
      for(size_t i=0 ; i<N ; ++i){
        tmp_v[i] = data[j][index[i]];
      }
      swap(data[j],tmp_v);
    }
  }

  size_t number_of_rows(){return data[0].size();}
  size_t number_of_columns(){return data.size();}
  
  std::vector<std::vector<T> > data;
private:
  std::map<std::string,int> datamap;
  std::vector<std::string> column_names;
  std::string filename;
};

// this is a summation algorithm that maintains percision for large sequences of numbers
template <typename T>
double KleinSum(std::vector<T> &input){
  double s = 0.0,cs = 0.0,ccs = 0.0,c,cc;
  size_t n = input.size();
  for(size_t i = 0 ; i<n ; ++i){
    double t = s + input[i];
    if( abs(s) >= abs(input[i]) ){
      c = (s - t) + input[i];
    }else{
      c = (input[i] - t) + s;
    }
    s = t;
    t = cs + c;
    if( abs(cs) >= abs(c) ){
      cc = (cs - t) + c;
    }else{
      cc = (c - t) + cs;
    }
    cs = t;
    ccs = ccs + cc;
  }
  return s + cs + ccs;
}


template <typename T,typename F>
double PairWiseSum(T *begin,T *end,F value){
  double sum = 0;
  size_t n = end - begin;
  if(n <= 10){
    for(T* p=begin ; p != end ; ++p){
      sum += value(p);
    }
  }else{
    sum = PairWiseSum<T,F>(begin,begin + n/2,value) + PairWiseSum<T,F>(begin + n/2,end,value);
  }
  return sum;
}

/// Does a parwise sumation of a vector which incresses the precision of the sumation to ~O(epsilon log(n) )
template <typename T>
double PairWiseSum(T *begin,T *end){
  double sum = 0;
  size_t n = end - begin;
  if(n <= 10){
    for(T* p=begin ; p != end ; ++p){
      sum += *p;
    }
  }else{
    sum = PairWiseSum(begin,begin + n/2) + PairWiseSum(begin + n/2,end);
  }
  return sum;
}

/// Does a parwise sumation of a vector which incresses the precision of the sumation to ~O(epsilon log(n) )
template <typename T>
double PairWiseSum(std::vector<T> &v){
  return PairWiseSum(v.data(),v.data() + v.size());
}

/// This version allows you to specify a function( value(T * p) ), that returns the value of a pointer to a T type that is to be summed.  This is useful for summing the squares or summing some a particular variable within a list of objects.  A lambda function is particularly useful here.
template <typename T,typename F>
double PairWiseSum(std::vector<T> &v,F value){
  return PairWiseSum(v.data(),v.data() + v.size(),value);
}

/// class for adding large amounts of numbers with less error than the simple sum
template <typename T>
class SUMMER{

public:

  SUMMER():batch(100),n(0),current(0.0),ntotal(0){};
  SUMMER(size_t batchsize)
  :batch(batchsize),n(0),current(0.0),ntotal(0){};

  /// add another number
  void operator+=(T x){
    ++n;
    current += x;
    ++ntotal;
    if(n % batch == 0){
      v.push_back(current);
      n=0;
      current = 0;
    }
  }
  
  /// returns the current total
  T operator*(){
    if(v.size() ==0 ) return 0;
    v.push_back(current);
    return Utilities::PairWiseSum(v);
  }
  
  /// reset to start over, frees memory
  void reset(){
    std::vector<T> dump;
    std::swap(v,dump);
    n = ntotal = 0;
    current = 0;
  }
  
  /// returns the number of numbers that have been added
  size_t total_number(){
    return ntotal;
  }
  
private:
  size_t batch;
  size_t n;
  size_t ntotal;
  T current;
  std::vector<T> v;
};

}  // Utilities

#endif
