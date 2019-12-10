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


#include "core/point.hxx"
#include "core/options.hxx"
#include <math.h>
#include <stdio.h>


Pool* Point :: pool = new Pool("Point",sizeof(Point),128) ;
// Creates the class variable 'pool', where all unused points will be
// stored, waiting for reuse (pseudo-dynamic allocation).


Point :: Point (real X, real Y, real Z)
    : GarbagedObject ()
{
    // Constructor. Initializes the receiver to a point with coordinates (X,Y,Z).

    x = X ;
    y = Y ;
    z = Z ;
}


Point* Point :: Add (Point* q)
{
    // p = p + q. Returns p (the receiver).

    x = x + q->x ;
    y = y + q->y ;
    z = z + q->z ;

    return this ;
}


real Point :: GetCoordinate (int i)
{
    // Returns the i-th coordinate of the receiver.

# ifdef REQUIRE
    Require("'i' is in [1..3]", i >= 1 && i <= 3) ;
    InMethod("Point::GetCoordinate(i)") ;
# endif

    if (i == 1)
	return x ;
    else if (i == 2)
	return y ;
    else
	return z ;
}


real Point :: GetDistanceTo (Point* p)
{
    // Returns the distance between the receiver and 'p'.

    Point *temp ;
    real  answer ;

    temp   = p->Minus(this) ;
    answer = temp->GetNorm() ;
    delete temp ;
  
    return answer ;
}


real Point :: GetNorm ()
{
    // Returns the norm of the coordinate array of the receiver.

    return sqrt (x*x + y*y + z*z) ;
}


void Point :: Print ()
{
    // A debugging tool. Prints the receiver on standard output.
  
    printf ("Point (% 6.4f,% 6.4f,% 6.4f)\n",x,y,z) ;
}


boolean Point :: IsIdenticalTo (Point* p)
{
    // Returns True if the receiver and 'p' have identical coordinates (within
    // a certain tolerance).

    return (fabs(x - p->x) < DEFAULT_TOLERANCE &&
	    fabs(y - p->y) < DEFAULT_TOLERANCE &&
	    fabs(z - p->z) < DEFAULT_TOLERANCE);
}


Point* Point :: Normalized ()
{
    // Returns the receiver after division by its 2-norm.

    real norm ;

    norm = GetNorm() ;

# ifdef CHECK
    if (norm == ZERO)     // DEFAULT_TOLERANCE is too large
	Error("Point::Normalized","norm is 0") ; 
# endif

    return Multiply(ONE/norm) ;
}


Point* Point :: Subtract (Point* q)
{
    // p = p - q. Returns p (the receiver).

    x = x - q->x ;
    y = y - q->y ;
    z = z - q->z ;

    return this ;
}


void Point :: PrintIn(FILE *fp)
{
    fprintf(fp,"%lf %lf %lf\n",x,y,z);
}
