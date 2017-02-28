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
#include <deque>

#endif
#include <mutex>

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
      reference operator*() const { return (reference)(**it); }
      
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
      reference operator[](difference_type n) const { return (reference)*it[n]; }
      
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
      const reference operator*() const { return (reference)*it; }
      
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
  

  /// class for keeping track of the range of a variable
  template<typename T>
  class Range{
  public:
    Range(T init_max,T init_min){
      range[0] = std::max(init_min,init_max);
      range[1] = std::min(init_min,init_max);
    }
    Range(T first_value){
      range[0] = range[1] = first_value;
    }
    
    void update(const T &val){
      range[0] = range[0] < val ? range[0] : val;
      range[1] = range[1] > val ? range[1] : val;
    }
    
    T max(){return range[1];}
    T min(){return range[0];}
    
  private:
    T range[2];
  };
  
  
#ifdef ENABLE_CLANG
  /// This is a class for generating random numbers. It is actually just a wrapper for some std random classes.
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
    
  };
  
  /// Shuffles a vector into a random order
  template <typename T, typename R>
  void shuffle(
               std::vector<T> &vec   /// The vector to be shuffled
               ,R ran               /// a random number generator so that ran() gives a number between 0 and 1
  ){
    T tmp;
    size_t ran_t;
    if(vec.size() < 2) return;
    for (size_t i = vec.size()-1; i>0; --i) {
      ran_t = (size_t)(ran()*(i+1));
      tmp = vec[ran_t];
      vec[ran_t] = vec[i];
      vec[i] = tmp;
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
  
  
  /// This should reorder the vector in place according to the permutation index p
  template <typename T>
  void apply_permutation_in_place(
                                  std::vector<T>& vec,
                                  const std::vector<std::size_t>& p)
  {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
      if (done[i])
      {
        continue;
      }
      done[i] = true;
      std::size_t prev_j = i;
      std::size_t j = p[i];
      while (i != j)
      {
        std::swap(vec[prev_j], vec[j]);
        done[j] = true;
        prev_j = j;
        j = p[j];
      }
    }
  }
  /// This should reorder the vector in place according to the permutation index p
  template <typename T>
  void apply_permutation_in_place(
                                  std::deque<T>& vec,
                                  const std::vector<std::size_t>& p)
  {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
      if (done[i])
      {
        continue;
      }
      done[i] = true;
      std::size_t prev_j = i;
      std::size_t j = p[i];
      while (i != j)
      {
        std::swap(vec[prev_j], vec[j]);
        done[j] = true;
        prev_j = j;
        j = p[j];
      }
    }
  }
  
#ifdef ENABLE_FFTW
  /** \brief Calculates power spectrum from a 2d map or the cross-power spectrum between two 2d maps.
   *
   *   Adaptied from Carlo Gioccoli's pl() routine.
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
                       ,int zeropaddingfactor
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
  
  /** \brief Read in data from an ASCII file with two columns
   */
  template <class T1,class T2>
  void read2columnfile(
                       std::string filename    /// input file name
                       ,std::vector<T1> &x     /// vector that will contain the first column
                       ,std::vector<T2> &y     /// vector that will contain the first column
                       ,std::string delineator = " "  /// specific string the seporates columns, ex. ",", "|", etc.
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
                       ,std::vector<T2> &y     /// vector that will contain the first column
                       ,std::vector<T3> &z     /// vector that will contain the first column
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
  
  
  /** \brief This is a thread-safe container wrapper.
   A reference to it can be passed to multiple threads without
   worrying about them interfering.  It is mostly meant as dynamical
   push only storage with STL data conatiners.  The functions mirror those
   of the STL data conatiners.
   */
  template<typename Container>
  class LockableContainer
  {
  public:
    typedef typename Container::value_type value_type;
    //    typedef typename Container::reference_type reference_type;
    typedef typename Container::value_type& reference_type;
    
    LockableContainer(){}
    ~LockableContainer(){}
    
    // make it uncopyable
    LockableContainer(const LockableContainer&) = delete;
    LockableContainer& operator=(const LockableContainer&) = delete;
    
    void push_back(value_type &t) {
      std::lock_guard<std::mutex> lock(m);
      container.push_back(t);
    }
    
    void push_front(value_type &t) {
      std::lock_guard<std::mutex> lock(m);
      container.push_front(t);
    }
    
    void resize(size_t N){
      std::lock_guard<std::mutex> lock(m);
      container.resize(N);
    }
    
    void reserve(size_t N){
      std::lock_guard<std::mutex> lock(m);
      container.reserve(N);
    }
    
    void clear(){
      std::lock_guard<std::mutex> lock(m);
      container.clear();
    }
    
    value_type operator[](size_t i) const{
      return container[i];
    }
    
    size_t capacity() const{
      return container.capacity();
    }
    size_t size() const{
      return container.size();
    }
    
    /// swap contents of conatiners
    void swap(Container &output){
      std::lock_guard<std::mutex> lock(m);
      std::swap(output,container);
    }
    
  private:
    Container container;
    std::mutex m;
  };
  
  /** \brief Class that make constructing and using a linear numerical lookput table easier.
   
   Uses linear interpolation.  std::bind() can be useful to get a function or method into a form that
   will be accepted by the constructor.
   **/
  template<class T>
  struct LinearLookUpTable{
  public:
    LinearLookUpTable(std::function<T(T)> func    /// function to sample from
                      ,T my_xmin                     /// minimum x value of table
                      ,T my_xmax                     /// maximum x value of table
                      ,size_t N                   /// number of points in table
                      ):xmax(my_xmax),xmin(my_xmin)
    {
      Utilities::fill_linear(x,N,xmin,xmax);
      y.resize(N);
      for(size_t i = 0;i<N; ++i){
        y[i] = func(x[i]);
      }
      dx = x[1] - x[0];
    }
    
    /// get linear interpolate at my_x, returns false if out of bounds
    bool operator()(T my_x,T &my_y){
      
      if(my_x < xmin || my_x >= xmax){   // out of bounds
        if(xmax == my_x){
          my_y = y.back();
          return true;
        }
        return false;
      }
      
      size_t i = (size_t)( (my_x - xmin)/dx );
      my_y = y[i] + (y[i+1] - y[i])*(my_x - x[i])/(x[i+1] - x[i]);
      return true;
    }
    
  private:
    std::vector<T> y;
    std::vector<T> x;
    double dx;
    double xmin;
    double xmax;
  };
  
  /** \brief Class that make constructing and using a linear numerical lookput table easier.
   Uses linear interpolation on the log values of both x and y.
   
   std::bind() can be useful to get a function or method into a form that
   will be accepted by the constructor.
   **/
  template <class T>
  class LogLookUpTable{
  public:
    LogLookUpTable(std::function<T(T)> func    /// function to sample from
                   ,T xmin                     /// minimum x value of table
                   ,T xmax                     /// maximum x value of table
                   ,size_t N                   /// number of points in table
    ){
      Utilities::fill_linear(x,N,log(xmin),log(xmax));
      y.resize(N);
      for(size_t i = 0;i<N; ++i){
        y[i] = log(func( exp(x[i]) ));
      }
      dx = x[1] - x[0];
      xmin = x[0];
      xmax = x.back();
    }
    
    bool operator()(T my_x,T &my_y){
      T lgx = log(my_x);
      
      if(lgx < xmin || lgx >= xmax){   // out of bounds
        if(xmax == lgx){
          my_y = exp( y.back() );
          return true;
        }

        return false;
      }
      
      size_t i = (size_t)( (lgx - xmin)/dx );
      
      my_y = exp( y[i] + (y[i+1] - y[i])*(lgx - x[i])/(x[i+1] - x[i]) );
      return true;
    }
    
  private:
    std::vector<T> y;
    std::vector<T> x;
    double dx;
    double xmin;
    double xmax;
  };
  
  /// cubic b-spline kernel in different dimensions normalized to 1
  template <int d>
  double Bspline(double q){
  }
  template <>
  inline double Bspline<1>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return 1.5*(0.25*q2*q2*q2 - q1*q1*q1);
    return q2*q2*q2/6;
  }
  template <>
  inline double Bspline<2>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return (0.25*q2*q2*q2 - q1*q1*q1)*0.3183098861837907;
    return q2*q2*q2*0.07957747154594767;
  }
  template <>
  inline double Bspline<3>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return (0.25*q2*q2*q2 - q1*q1*q1)/pi;
    return q2*q2*q2*0.07957747154594767;
  }
}
#endif
