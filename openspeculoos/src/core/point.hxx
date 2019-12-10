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

#ifndef POINT
#define POINT

/*!
 * \file point.hxx
 * \brief Manages point types
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/garbaged.hxx"
#include "core/pool.hxx"
#include "core/util.hxx"

/*! \class Point
   * \brief Manages 3D points.
   * 
   * A point is described by its 3 coordinates <{x}>, <{y}>, <{z}>.
   * In 2D problems the z-ccordinate is set to 0.
   * 
   * When the user creates a geometrical mesh, he/she usually starts by
   * creating some points, then vertices based on each point, then edges
   * based on these vertices, then faces based on these edges, etc.
   * 
   * Profiling tests showed that points are very frequently created and
   * deleted, so a pool of points (the class attribute <{pool}>) is used to
   * store available points.
   */ 

class Point : public GarbagedObject
{

protected :

    real x ;
    real y ;
    real z ;

    static Pool* pool ;

public :

    // = PSEUDO-DYNAMIC ALLOCATION

    inline void*   operator new (size_t) ;
    inline void    operator delete (void* p) ;

    // = CONSTRUCTORS

    Point (real X=ZERO, real Y=ZERO, real Z=ZERO) ;   // should be in-lined

    // = COORDINATES

    real           GetCoordinate (int i) ;
    inline real    GetX () ;
    inline real    GetY () ;
    inline real    GetZ () ;

    // = OPERATORS (OBSOLETE)
 
    inline Point   operator -  (Point& p) ;
    inline Point   operator %  (Point& p) ;

    // = OPERATIONS

    Point*  Add (Point* p) ;
    inline Point*  Minus (Point* p) ;
    inline Point*  Plus (Point* p) ;
    inline Point*  Multiply (real alpha) ;
    Point*  Subtract (Point* p) ;
    inline Point*  Times (real alpha) ;
    inline real    Dot (Point* p) ;

    // = COMPARISON

    boolean        IsIdenticalTo (Point* p) ;

    // = COPY

    inline Point*  Duplicate () ;

    // = NORM 

    real           GetNorm () ;
    Point*         Normalized () ;
    real           GetDistanceTo (Point* p) ;
 
    // = PRINTING

    void           PrintIn(FILE *fp);
    void           Print () ;
} ;


//____________________________ Point (inline) _________________________________


void* Point :: operator new (size_t)
{
    // Pseudo-dynamic allocation of a point.

    return pool->Alloc() ;
}


void Point :: operator delete (void* p)
{
    // Pseudo-dynamic allocation of 'p'.

    pool->Free(p) ;
}


Point Point :: operator - (Point& p) 
{
    // Returns a point whose coordinates are those of the receiver minus those of
    // 'p'.
 
    return Point (x-p.x, y-p.y, z-p.z) ;
}


Point Point :: operator % (Point& p)
{
    // Exterior product of the receiver and 'p'.

    return Point(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x) ; 
}


real Point :: Dot (Point* p) 
{
    // Returns the scalar product of the receiver and 'p'.
 
    return x*p->x + y*p->y + z*p->z ; 
}   


Point* Point :: Duplicate ()
{
    // Returns a copy of the receiver.

    return new Point(x,y,z) ;
}


real Point :: GetX () 
{
    // Returns the first coordinate of the receiver.

    return x ;
}


real Point :: GetY () 
{
    // Returns the second coordinate of the receiver.

    return y ;
}


real Point :: GetZ () 
{
    // Returns the third coordinate of the receiver.

    return z ;
}


Point* Point :: Minus (Point* p)
{
    // r = p - q. Returns r (a new Point)

    return Duplicate()->Subtract(p) ;
}


Point* Point :: Plus (Point* p)
{
    // r = p + q. Returns r (a new Point)

    return Duplicate()->Add(p) ;
}


Point* Point :: Multiply (real a)
{
    // p = p * a. Returns p (the receiver).

    x = a*x ;
    y = a*y ;
    z = a*z ;

    return this ;
}


Point* Point :: Times (real a)
{
    // q = p * a. Returns q (a new point).

    return new Point (x*a,y*a,z*a) ;
}


#endif
