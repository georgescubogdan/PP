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

#ifndef POOL
#define POOL

/*!
 * \file pool.hxx
 * \brief Manages pseudo-dynamic allocation of objects
 * \author {Lots of authors}
 * \version 1.0
 */

#include <stdio.h>
#include <string>

using namespace std ;

/*! \class Pool
   * \brief Manages 3D points.
   * 
   * A pool is a list of objects of the same type (say, of type Point)
   * which are available for use. The points can be requested from, and
   * released to, the pool repeatedly, so repeated dynamic allocation/
   * deallocation are avoided.
   * 
   * In order to create/delete objects, calls to the default C++ new/delete
   * operators (which dynamically allocate/deallocate memory through calls
   * to the operating system) are replaced by calls to overloaded new/
   * delete operators (which call methods Alloc/Free of class Pool, thus
   * avoiding dynamic allocation/deallocation) of class Point.
   * As far as the user is concerned, the pool is transparent: a point is
   * "created" by means of
   *    point = new Point() ;
   * and "deleted" by means of
   *    delete point ;
   * as if there where no pool.
   * 
   * Note that allocating objects from a pool avoid allocating them 
   * dynamically, but it does not avoid invoking the constructor.
   * 
   * At execution end, the pool can be deleted. However, since the pool
   * takes no record of the object it allocated (in method <{Grow}>), it
   * cannot dynamically deallocate the objects. This memory leak is usually
   * not detected by Purify. Therefore, in order to check that the user
   * does not forget to release the objects he/she obtained from the pool,
   * the pool is provided with an additional attribute <{used}>, which 
   * gives at any time the number of objects currently used by the user.
   * When the pool is deleted, this number should be zero (and the pool
   * should thus be full).
   * 
   * Attribute <{name}> contains typically the type of the objects, e.g.,
   * "Point".
   * 
   * The <{Link}> structure uses the memory occupied by the object to store
   * its attribute <{next}>.
   * Attribute <{nelem}> contains the number of objects (eg, of points) to
   * be allocated every time the pool is empty (default: <{nelem}> = 1).
   * 
   * See Stroustrup, The C++ programming language (2nd ed.), Addison-
   * Wesley, 1991, pp 473-474 (and also pp 176-177).
   */ 
class Pool
{

protected :

    struct Link { Link* next ; } ;

    string name ;
    int    esize ;
    Link*  head ;
    int    nelem ;
    int    used ;
  
    void Grow() ;

public :

    Pool (string aname, int sz, int n) ;
    ~Pool () ;

    inline void* Alloc () ;
    inline void  Free (void* obj) ;
} ;


//______________________________ Pool (inline) ________________________________


void* Pool :: Alloc ()
{
    // Returns a pointer to an object. Enlarges the receiver if necessary.

    Link *p ;

    if (head == NULL)
	Grow() ;

    p    = head ;
    head = p->next ;

    used++ ;

    return p ;
}


void Pool :: Free (void* obj)
{
    // Dumps 'obj' into the receiver. 'obj' thus becomes available for future
    // reuse.

    Link *p ;

    p       = (Link*)obj ;
    p->next = head ;
    head    = p ;

    used-- ;
}


#endif
