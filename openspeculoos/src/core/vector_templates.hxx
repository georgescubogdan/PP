/*
* 
* This file is part of OpenSPECULOOS
* Copyright (c) 2008 Ecole Polytechnique Federale de Lausanne (EPFL)
*
* Contact: Nicolas Fietier, EPFL, Station 9, 1015 Lausanne, Switzerland
* E-mail:  nicolas.fietier@epfl.ch
*
* OpenSPECULOOS is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* OpenSPECULOOS is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with OpenSPECULOOS.  If not, see <http://www.gnu.org/licenses/>.
* 
*/

#ifndef VECTOR_TEMPLATES
#define VECTOR_TEMPLATES

/*!
 * \file templates.hxx
 * \brief Templates
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/util.hxx"
#include "core/options.hxx"
#include <stdio.h>
#include <math.h>
#include <string>
#include <stdlib.h>

/*! \brief Manages vectors of any type.
   * 
   * 
   * Implements an array containing values of type T. T may be anything:
   * basic-type (like 'int'), pointer to basic-type (like 'int*'), pointers
   * to objects (like 'Element*'), etc.
   * 
   * Possibility to resize dynamically the array is provided. Search
   * capability is also available together with simple affectations and
   * equality tests.
   * 
   * Since resizing dynamically an array may be cost-ineffective, an array
   * is supplied with two size attributes:
   *   . <{allocatedSize}> is the physical size of the vector, i.e., the
   *     number of coefficients of the array <{values}>;
   *   . <{size}> is the logical size of the vector. At any time, <{size}>
   *     can be modified to any value smaller than or equal to 
   *     <{allocatedSize}>, without performing dynamic reallocation.
   * 
   * Bounds are checked if the REQUIRE compilation option is defined in
   * file "options.hxx".
   * 
   * The default value of <{isAlias}> is 'false'. A vector v2 can be 
   * aliased to another vector v1. This means that attribute <{values}> of
   * v2 is set to attribute <{values}> of v1, so both vectors share the 
   * same values. Aliasing is useful for handling a one-dimensional vector
   * array as a two-dimensional matrix array.
   * If v2 is deleted, in order to prevent it from deleting its <{values}>,
   * its flag <{isAlias}> is set to 'true'.
   */ 
template <class T> class Vector
{ 

protected :

    int     allocatedSize ;
    int     size ;
    T*      values ;
    boolean isAlias ;

    void    Alloc () ;
    boolean ValidIndex (int i) ;

public:

    // = CONSTRUCTORS

    Vector () ;
    Vector (int n) ; 
    ~Vector () ;

    // = INQUIRY

    inline T&       At (int i) ;
    int             GetIndexOf (T val) ;
    inline int      GetSize () ;
    inline boolean  Has (T val) ;
    inline boolean  HasIndex (int i) ;
    boolean         HasOnly (T val) ;
    boolean         IsIdenticalTo (Vector<T>* x) ;

    // = COPY

    void            CopyFrom (Vector<T>* x) ;

    // = OPERATIONS

    T*              GetValues () ;
    void            DefineAlias (Vector<T>* v, int start, int length) ;
    Vector<T>*      Duplicate () ;
    void            Put (T val) ;
    int             PutAndGetIndex(T val);

    void            Remove (T val) ;
    void            RemoveIndex (int i) ;
    void            Resize (int n) ;
    void            SetValues (T val) ;
    inline void     SetSize (int n) ;
} ;


class  Communicator;
class  DiagonalMatrix;
class  DirichletCondition ;
class  Distribution ;
class  Edge ;
class  Element ;
class  ElementaryField;
class  Field;
class  FlatField ;
class  FullMatrixLU ;
class  Matrix ;
class  TensorMatrix ;
class  NeumannCondition ;
class  OperatorQ1D ;
class  ParentElement ;
class  ParentElementGL ;
class  ParentElementGLL ;
class  ParentElementGLLI ;
class  ParentElementHP ;
class  ParentElementIGLL ;
class  Quad;
class  Vertex;
class  RealVector;
class  Point;
class  Mesh;
class  Mesh2D;
class  Mesh1D;

template class Vector<Point*> ; 
template class Vector<Mesh2D*>;
template class Vector<Mesh*>;
template class Vector<Mesh1D*>;

template class Vector<Communicator*> ;
template class Vector<DiagonalMatrix*> ;  
template class Vector<DirichletCondition*> ;
template class Vector<Distribution*> ;
template class Vector<Edge*> ;
template class Vector<Element*> ;
template class Vector<Element**> ;
template class Vector<ElementaryField*>;
template class Vector<FlatField*> ;
template class Vector<FullMatrixLU*> ;
template class Vector<Field*> ;
template class Vector<Matrix*> ;
template class Vector<NeumannCondition*> ;
template class Vector<OperatorQ1D*> ;
template class Vector<ParentElement*> ;
template class Vector<ParentElementGL*> ;
template class Vector<ParentElementGLL*> ;
template class Vector<ParentElementGLLI*> ;
template class Vector<ParentElementHP*> ;
template class Vector<ParentElementIGLL*> ;
template class Vector<Quad*> ;
template class Vector<RealVector*> ;
template class Vector<TensorMatrix*> ;
template class Vector<Vertex*> ;
template class Vector<boolean> ;
template class Vector<int *>;
template class Vector<real>;
template class Vector<string>;

//____________________________ Vector (inline) ________________________________


template <class T>  
T& Vector<T> :: At (int i) 
{
    // Returns the 'i'-th component of the receiver.

# ifdef REQUIRE 
    Require("valid index", ValidIndex(i)) ;
    InMethod("Vector::At(i)") ;
# endif

    return values[i-1] ; 
}


template <class T> 
T* Vector<T> :: GetValues ()
{
    // Returns a pointer to the first value of the receiver.
    // This method is an encapsulation crime: avoid using it, as much as 
    // possible!

    return values ;
}


template <class T> 
int Vector<T> :: GetSize () 
{
    // Returns the size of the receiver.

    return size ; 
}


template <class T> 
boolean Vector<T> :: Has (T val)
{
    // Returns True if the receiver contains 'val', else returns False.
    // If T is a pointer type, only pointer equality is tested, not object
    // equality.

    return (GetIndexOf(val) != 0) ;
}


template <class T> 
boolean Vector<T> :: HasIndex (int i)
{
    // Returns True if 'i' is within the receiver's index range.

    return (i >= 1 && i <= size) ;
}


template <class T>
void Vector<T> :: SetSize (int newSize)
{
    // Resets the receiver's size to 'newSize'. The values are not affected.
    // The receiver is not allowed to grow.

# ifdef REQUIRE
    Require("'newSize' is positive", newSize >= 0) ;
    Require("'newSize' not too large", newSize <= allocatedSize) ;
    Require("vector is not aliased", ! isAlias) ;
    InMethod("Vector::SetSize(newSize)") ;
# endif

    size = newSize ;
}

template <class T>
Vector<T> :: Vector ()
{
    // Constructor. Initializes the receiver to a vector of size 0.

    allocatedSize = 0 ;
    size          = 0 ;
    isAlias       = false ;

    Alloc() ;
}


template <class T>
Vector<T> :: Vector (int n)
{
    // Constructor. Initializes the receiver to a vector of size 'n'.
    // The values are not initialized (contain garbage).

# ifdef REQUIRE
    Require("n >= 0", n >= 0) ;
    InMethod("Vector(n)") ;
# endif

    allocatedSize = n ;
    size          = n ;
    isAlias       = false ;

    Alloc() ;
}


template <class T>
Vector<T> :: ~Vector ()
{
    // Destructor.
    // If the receiver is an alias of another vector, does nothing.
    // If the receiver is a standard vector, deletes the array; the values are
    // not deleted (if of pointer type).
    //
    // Programming note: do not attempt to define an additional attribute
    //    'deletableValues' which, when set to True, would also delete the
    //    values within 'values'. The reason is that the 'delete values[i]'
    //    instruction is not valid (error message at link time), because the
    //    class is a template<T>, not a template<T*>. 
    //    If the values are pointers to object, you can delete them, but not
    //    from this class, only from your class.

    if (! isAlias)
	delete [] values ;
}


template <class T> 
void Vector<T> :: Alloc ()
{
    // Allocates storage for the receiver.

# ifdef REQUIRE
    Require("allocated size >= 0", allocatedSize >= 0) ;
    InMethod("Vector::Alloc()") ;
# endif

    if (allocatedSize > 0) {
	values = new T[allocatedSize] ;
#   ifdef CHECK
	if (! values)
	    Error("Vector::Alloc()","free store exhausted") ;
#   endif
    }
    else
	values = NULL ;
}


template <class T> 
void Vector<T> :: CopyFrom (Vector<T>* v)
{
    // Copies 'v' into the receiver. It is a shallow copy, i.e., if the values
    // are pointer-type values, the pointed objects are not copied.

# ifdef REQUIRE
    Require ("same size", v->size == size) ;
    InMethod("Vector::CopyFrom(v)") ;
# endif

    int i ;

    for (i=0 ; i<size ; i++)
	values[i] = v->values[i]  ;
}


template <class T> 
void Vector<T> :: DefineAlias (Vector<T>* v, int start, int length)
{
    // Defines the receiver as an alias to the part of 'v' starting from index
    // 'start' and having length 'length'.

# ifdef REQUIRE
    Require("valid argument 'start'", v->ValidIndex(start)) ;
    Require("valid argument 'length'",length>0 && v->ValidIndex(start+length-1));
    InMethod("Vector::DefineAlias(v,start,length)") ;
# endif

    if (! isAlias)
	delete [] values ;
    isAlias = true ;

    values = v->values + start - 1 ;
    size   = length ;

    // allocatedSize is not defined
}


template <class T> 
Vector<T>* Vector<T> :: Duplicate ()
{
    // Returns a copy of the receiver. It is a shallow copy, i.e., if the values
    // are pointer-type values, the pointed objects are not duplicated.

    Vector<T>* answer ;
    int        i ;

    answer = new Vector<T>(size) ;
    for (i=0 ; i<size ; i++)
	answer->values[i] = values[i] ;

    return answer ;
}


template <class T>
int Vector<T> :: GetIndexOf (T val)
{
    // Returns the index (between 1 and 'size') of the first occurence of 'val' 
    // in the receiver. Returns 0 if cannot find 'val'.

    int i ;

    for (i=0 ; i<size ; i++)
	if (values[i] == val)
	    return i+1 ;

    return 0 ;
}

template <class T>
boolean Vector<T> :: HasOnly (T val) 
{
    // Returns 'true' if all values of the receiver are 'val', else returns
    // 'false'.

    int i ;

    for (i=0 ; i<size ; i++)
	if (values[i] != val)
	    return false ;

    return true ;
}


template <class T>
boolean Vector<T> :: IsIdenticalTo (Vector<T>* v) 
{
    // Returns 'true' if the receiver and 'v' have the same contents.

    int i, answer ; 

    if (size != v->GetSize()) 
	answer = false ;
    else {
	answer = true ;
	for (i=0 ; i<size ; i++) 
	    if (values[i] != v->values[i]) {
		answer = false ;
		break ;
	    }
    }

    return answer ;
}


template <class T>
void Vector<T> :: Put (T val)
{
    // Appends 'val' to the receiver, i.e., sets the (size+1)-th entry of the
    // receiver to 'val'. Reallocates the array of values if necessary.

    Resize(size+1) ;
    values[size-1] = val ;
}

template <class T>
int Vector<T> :: PutAndGetIndex(T val)
{
    // Appends 'val' to the receiver, i.e., sets the (size+1)-th entry of the
    // receiver to 'val'. Reallocates the array of values if necessary.

    Resize(size+1) ;
    values[size-1] = val ;

    return size;
}

template <class T>
void Vector<T> :: Remove (T x)
{
    // Removes all occurences of 'x' in the receiver. If found 'x', compresses
    // the receiver; the size of the receiver is decremented by the number of
    // occurences of 'x'. The array of values is not reallocated.
 
    int i ;

    while ((i=GetIndexOf(x)) != 0)
	RemoveIndex(i) ;
}


template <class T>
void Vector<T> :: RemoveIndex (int i)
{
    // Removes the i-th component of the receiver. Compresses the receiver; the
    // size of the receiver is decremented by 1. The array of values is not 
    // reallocated.

# ifdef REQUIRE
    Require("valid index", ValidIndex(i)) ; 
    InMethod("Vector::RemoveIndex(i)") ;
# endif

    int j ;

    for (j=i ; j<allocatedSize ; j++) 
	values[j-1] = values[j] ;

    size-- ;
}


template <class T>
void Vector<T> :: Resize (int newSize)
{
    // Sets to 'newSize' the size of the receiver.
    //
    // If 'newSize' is larger than the current size, the new values (from index 
    // "size+1" to index "newSize") are initialized to 0 (or NULL).
    //
    // Reallocates the array of values if 'newSize' is larger than the allocated
    // size.

# ifdef REQUIRE
    Require("'newSize' is positive", newSize > 0) ;
    Require("vector is not aliased", ! isAlias) ;
    InMethod("Vector::Resize(newSize)") ;
# endif

    T         *oldValues ;
    int       i, oldSize ;
    const int SIZE_INCREMENT = 16 ;

    oldSize = size ;
    if (newSize > allocatedSize) {          // need to reallocate 
	oldValues     = values ;
	size          = newSize ;
	allocatedSize = max(newSize,allocatedSize+SIZE_INCREMENT) ;
	Alloc() ;
	for (i=0 ; i<oldSize ; i++)
	    values[i] = oldValues[i] ;
	delete [] oldValues ;
    }
    else                                    // no need to reallocate
	size = newSize ;

    for (i=oldSize ; i<newSize ; i++)
      values[i] = T() ;
}


template <class T>
void Vector<T> :: SetValues (T val)
{
    // Sets every component of the receiver to 'val'.

    int i ;

    for (i=0 ; i<size ; i++)
	values[i] = val ;
}


template <class T>
boolean Vector<T> :: ValidIndex (int i)
{
    // Returns 'true' if 'i' is within {1,size}, else returns 'false'.

    if (i<1 || i>size) {
	printf ("\nWarning: using index %d in Vector of size %d\n",i,size) ;
	return false ;
    }
    else
	return true ;
}


#endif
