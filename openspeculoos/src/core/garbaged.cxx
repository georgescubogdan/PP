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

// garbaged.cxx

#include "core/garbaged.hxx"
#include <stdio.h>


void GarbagedObject :: BreakCyclicReferences ()
{
    // Default implementation: does nothing.
}


void GarbagedObject :: PrintCounter ()
{
    // Prints on standard output the number of objects ponting to the receiver.

    printf("  counter of object: %d\n",counter) ;
}


void ref (GarbagedObject* obj)
{
    // Increments by 1 the number of references to 'obj'.
    // Can be safely invoked if 'obj' is NULL.

    if (obj)
	(obj->counter)++ ;
}


void unref (GarbagedObject* obj)
{
    // Decrements by 1 the number of references to 'obj'. Deletes 'obj if it is
    // not referenced any longer.
    // Can be safely invoked if 'obj' is NULL.

    if (obj) {
	(obj->counter)-- ;
	obj->BreakCyclicReferences() ;
	if (obj->counter == 0)
	    delete obj ;
    }

}


