// File: List.H  List class for any data type.
#ifndef _List_H_
#define _List_H_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/array.h"

/*! \class List list.h oml/list.h
  \brief Numerical C-array like container with auto resizing.

  Lists have element indexes ranging from 0...n-1 and are indexed using the \c [i]
  syntax, just like C arrays.  The base index 0 cannot be changed. Lists support
  copy on write (COW) symmantics, which avoids deep copying (element by element)
  when they are returned by value from a function. One can add elements to the end of a List
  and it will automaticllay resize itself. Partial emulation of an STL \c std::vector<T>.

  \b Math:

  The special operators for List are (same as Array).
  - \c L*L does a direct multiply (not a dot product).
  - use \c Dot(L,L) to get a dot product.
  - \c L/L does a direct divide.

  \nosubgrouping
*/
template <class T> class List
  : public Indexable<T,List<T>,Full,Real,ArrayShape>
  , public Iterable<T,List<T> >
  , public TStreamableObject<List<T> >
{
  static const index_t start_size=16;
 public:
  /*! \name Constructors/Assignment
  Copy constructor and op=(List) are automatically supplied by the compiler.

  */
  //@{

  //! Empty list.
           List(                          ) : itsData(start_size), itsNext(0   ) {};
  //! Unitialized list.
  explicit List(index_t size              ) : itsData(size      ), itsNext(size) {};
  //! List filled with one value.
  List(index_t size,const T& fill) : itsData(size      ), itsNext(size) {Fill(*this,fill);}
  //! Allows construction from an expression template.
  template <class T1,class A> List(const Indexable<T1,A,Full,Real,ArrayShape>&);
  //! Allows assignment from a list.
  List& operator =(const List&);
  //! Allows assignment from an expression template.
  template <class T1,class A> List& operator =(const Indexable<T1,A,Full,Real,ArrayShape>&);
  //@}

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const List& l)
  {
    return os << static_cast<const TStreamableObject<List<T> >& >(l);
  }


  /*! \name Subscripting operators
    If \c DEBUG is defined, every index will be checked that it is in range.
    For fast write access with \c op[] make a \c Subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element access operator, fast and \e cannot trigger a COW operation.
  const T& operator[](index_t i) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
        T& operator[](index_t i)      ;
  //@}

  //! Returns number elements in the array. Should be the same as the allocated space.
  index_t size() const;
  //! Change the size of the array and optionally preserve as many values as possible.
  void    SetSize(index_t, bool preserve=false)      ;
  //! Does A[i]=A[index[i]] for i=0...n-1.  Used for sorting.
  void    ReIndex(const Array<index_t>& index);
  //! Extract a subarray.  Involves a deep copy.
  List    SubList(index_t, index_t) const;


  /*! \name List operatorations
    These functions try to emulate an STL std::vector.
  */
  //@{

  //! Add one element to the end of the List.
  void     push_back(const T&         )      ;
  //! Add one element to the end of the List.
  void     pop_back (                 )      ;
  //! Get last element.
  const T& back     (                 ) const;
  //! Get first element.
  const T& front    (                 ) const;
  //! Get last element.
        T& back     (                 )      ;
  //! Get first element.
        T& front    (                 )      ;
  //! Find first occurence of a value.
  index_t  find     (const T&         ) const;
  //! Insert one element in the middle. Expensive.
  void     InsertAt (const T&, index_t)      ;
  //! Remove one element in the middle. Expensive.
  void     RemoveAt (          index_t)      ;
  //! Locate value and remove the first occurence.
  bool     Remove   (const T&         )      ;
  //! Remove all occurences of a value.
  void     RemoveAll(const T&         )      ;
  //! Clear out the List. Size goes to zero.
  void     Empty    (                 )      ;
  //! Look if a value is in the list.
  bool     IsInList (const T&         ) const;

  List& operator&(const List& b);
  //@}

#if DEBUG
  #define CHECK(i) assert(i>=0&&i<itsSize)
#else
  #define CHECK(i)
#endif
  class ArraySubscriptor
  {
   public:
    ArraySubscriptor(Indexable<T,List,Full,Real,ArrayShape>& a)
      : itsPtr(static_cast<List&>(a).Get()), itsSize(a.size()) {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };
#undef CHECK
  /*! \name Subscriptors and Iterators.
    Iterators should be STL compatible.
   */
  //@{
  //! Read only iterator.
  typedef typename Iterable <T,List>::const_iterator  const_iterator ;
  //! Read/write iterator.
  typedef typename Iterable <T,List>::iterator iterator;
  //! Default Subscriptor.
  typedef ArraySubscriptor Subscriptor;
  //@}

 private:
  friend class ArraySubscriptor;
  friend class Iterable <T,List>;

  const T* Get() const;
        T* Get()      ;
  void  CheckIndex(index_t) const;
  void  SetCapacity(index_t, bool preserve=false); //Changes itsData but no itsNext;
  cow_array<T> itsData;
  index_t      itsNext;
};

//-------------------------------------------------------------------------------
//
//  These macros invoke List index bounds checking if DEBUG is on.
//
#if DEBUG
  #define CHECK(i) CheckIndex(i)
#else
  #define CHECK(i)
#endif

template <class T> inline const T&  List<T>::operator[](index_t i) const
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& List<T>::operator[](index_t i)
{
  CHECK(i);
  return itsData.Get()[i];
}

#undef CHECK

//-------------------------------------------------------------------------------
//
//  Some inline functions.
//
template <class T> inline index_t  List<T>::size() const
{
  return itsNext;
}

template <class T> inline const T* List<T>::Get() const
{
  return itsData.Get();
}

template <class T> inline T* List<T>::Get()
{
  return itsData.Get();
}

template <class T> inline bool  List<T>::IsInList(const T& t) const
{
  return find(t)<itsNext;
}

template <class T> template <class T1,class A> inline
List<T>::List(const Indexable<T1,A,Full,Real,ArrayShape>&x)
: itsData(x.size())
, itsNext(x.size())
{
  assert(itsNext <= itsData.size());
  ArrayAssign(*this,x);
}


template <class T> inline
List<T>& List<T>::operator=(const List<T>& x)
{
  itsData=x.itsData; //Shallow copy.
  itsNext=x.itsNext;
  assert(itsNext <= itsData.size());
  return *this;
}

template <class T> template <class T1,class A> inline
List<T>& List<T>::operator=(const Indexable<T1,A,Full,Real,ArrayShape>& x)
{
  SetSize(x.size());
  ArrayAssign(*this,x);
  assert(itsNext <= itsData.size());
  return *this;
}

template <class T> inline const T& List<T>::front() const
{
  return (*this)[0];
}

template <class T> inline const T& List<T>::back() const
{
  return (*this)[size()-1];
}

template <class T> inline T& List<T>::front()
{
  return (*this)[0];
}

template <class T> inline T& List<T>::back()
{
  return (*this)[size()-1];
}

template <class T> inline void List<T>::Empty()
{
  SetSize(start_size);
  itsNext=0;
}

//------------------------------------------------------------------------------#
//
//  Change the List to a new size.  As much old data as possible
//  is saved.
//
template <class T> void List<T>::SetSize(index_t newsize, bool preserve)
{
  SetCapacity(newsize,preserve);
  itsNext=newsize;
}

template <class T> void List<T>::SetCapacity(index_t newsize, bool preserve)
{
  assert(newsize >= 0);
  if (itsData.size()!=newsize)
  {
    index_t n=size();
    if (preserve && n>0)
    {
      List<T> dest(newsize);
      dest.itsNext=n;
      iterator d=dest.begin();
      for (const_iterator s=this->begin();s!=this->end();d++,s++) *d=*s;
      *this=dest;
    }
     else
    {
      *this=List<T>(newsize);
      itsNext=n;
    }
  }
}

template <class T> List<T> List<T>::SubList(index_t start, index_t stop) const
{
  assert(start>=0);
  assert(stop>=start);
  assert(stop<size());
  List<T> ret;
  for (index_t i=start;i<=stop;i++) ret.push_back((*this)[i]);
  return ret;
}

#if DEBUG
//------------------------------------------------------------------------------#
//
//  Check if an index is within bouds.
//
template <class T> void List<T>::CheckIndex(index_t i) const
{
  if (i>=itsNext || i<0)
  {
    OMLListIndexError(i,itsNext-1);
    assert(i<itsNext);
    assert(i>=0);
  }
}
#endif //DEBUG

template <class T> inline List<T> operator&(const List<T>& a, const List<T>& b)
{
  typedef typename List<T>::const_iterator CI;
  List<T> ret;
  for (CI ab=a.begin();ab!=a.end();ab++) ret.push_back(*ab);
  for (CI bb=b.begin();bb!=b.end();bb++) ret.push_back(*bb);
  return ret;
}

template <class T> inline void List<T>::ReIndex(const Array<index_t>& index)
{
  assert(size()==index.size());

  List<T> dest(size());

  Array<index_t>::const_iterator b=index.begin();
  iterator                       i=dest.begin();
  for (;b!=index.end();b++,i++) *i=(*this)[*b];
  *this=dest;
}


//---------------------------------------------------------------------------
//
//  Add one member to list and reallocate if nessecary.
//
template <class T> inline void List<T>::push_back(const T& element)
{
  index_t N=itsData.size();
  assert(itsNext <= N);
  if(itsNext == N) SetCapacity(N>0 ? N*2 : start_size,true);
  (*this)[itsNext++]=element;
}

//---------------------------------------------------------------------------
//
//  Pop one member off the back of the list.
//
template <class T> inline void List<T>::pop_back()
{
  assert(itsNext>0);
  itsNext--;
}

//---------------------------------------------------------------------------
//
//  Insert one object at specified location
//
template <class T> void List<T>::InsertAt(const T& element, index_t e)
{
  assert(e<itsNext);
  push_back(back());  //This takes care of any reallocation.
  Subscriptor s(*this);
  for (index_t i=size()-1;i>e; i--) s[i]=s[i-1];
  s[e]=element;
}

//---------------------------------------------------------------------------
//
//  Find and then delete an object
//
template <class T> inline bool List<T>::Remove(const T& element)
{
  index_t i=find(element);
  if (i<itsNext)
  {
    RemoveAt(i);
    return true;
  }
  return false;
}

//---------------------------------------------------------------------------
//
//  Remove all object that match the argument.
//
template <class T> void List<T>::RemoveAll(const T& element)
{
  while (Remove(element)) {};
}

//---------------------------------------------------------------------------
//
//  Remove one member of the list and dont bother to reallocate.
//
template <class T> void List<T>::RemoveAt(index_t e)
{
  assert(e<=size());
  Subscriptor s(*this);
  for (index_t i=e; i< size()-1; i++) s[i]=s[i+1];
  itsNext--;
}

//---------------------------------------------------------------------------
//
//  Search the list for an object.
//
template <class T> index_t List<T>::find(const T& element) const
{
  index_t i;
  for (i=0;i<size();i++) if ((*this)[i]==element) break;
  return i;
}

//! Append another list.
template <class T> inline List<T>& List<T>::operator&(const List<T>& b)
{
  typedef typename List<T>::const_iterator CI;
  for (CI bb=b.begin();bb!=b.end();bb++) this->push_back(*bb);
  return *this;
}

#endif //_List_H_
