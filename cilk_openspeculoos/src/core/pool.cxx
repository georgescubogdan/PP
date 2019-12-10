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

// pool.cxx

#include "core/pool.hxx"
#include "core/util.hxx"


Pool :: Pool (string aname, int sz, int n) 
{
    // Constructor. Initializes the receiver to a pool for objects of size 'sz'.
    // At every allocation, 'n' objects will be allocated at the same time.
    // 'sz' can be obtained by invoking sizeof(C), where C stands for the type of
    // the objects to be stored in the pool.

    name  = aname ;
    esize = (sz < sizeof(Link*) ? sizeof(Link*) : sz) ;
    head  = NULL ;
    nelem = n ;
    used  = 0 ;
}


Pool :: ~Pool ()
{
    // Destructor.
    // Unable to free the allocated chunks, because no record has been made of
    // the allocations in method Grow().

    if (used)
      printf("\nnumber of still used objects of Pool '%s': %d\n",name.c_str(),used) ;
}


void Pool :: Grow ()
{
    // Enlarges the receiver by allocating chunks of memory for new objects.
  
    char      *start, *last, *p ;
    const int chunk_size = esize ;

    start = new char[chunk_size*nelem] ;
    last  = &start[(nelem-1)*esize] ;
  
    for (p=start ; p<last ; p+=esize)
	((Link*)p)->next = (Link*)(p+esize) ;
    ((Link*)last)->next = NULL ;

    head = (Link*)start ;
}


