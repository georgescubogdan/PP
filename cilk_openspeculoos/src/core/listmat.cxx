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


#include "core/listmat.hxx"
#include "core/matrix.hxx"


ListMatrices :: ListMatrices ()
{
    // Constructor. Initializes the receiver to an empty list.

    parents  = new Vector<ParentElement*>() ;
    matrices = new Vector<Matrix*>() ;
}


ListMatrices :: ~ListMatrices ()
{
    // Destructor. The derivative matrices are deleted, the parent elements are
    // not.

    int i, n ;
  
    delete parents ;

    n = matrices->GetSize() ;
    for (i=1 ; i<=n ; i++)
	delete matrices->At(i) ;
    delete matrices ;
}


Matrix* ListMatrices :: Search (ParentElement* p)
{
    // Returns the matrix associated with 'p'. Returns NULL if there is no entry
    // for 'p'.

    int i, idx, n ;

    n = parents->GetSize() ;
    for (i=1 ; i<=n ; i++) {
	idx = parents->GetIndexOf(p) ;
	if (idx)
	    return matrices->At(idx) ;
    }

    return NULL ;
} 
  
  
void ListMatrices :: Store (ParentElement* p, Matrix* m)
{
    // Creates a new entry with 'p' and the associated matrix 'm'.
 
    int n ;

    n = parents->GetSize() ;
    parents->Resize(n+1) ;
    parents->At(n+1) = p ;
    matrices->Resize(n+1) ;
    matrices->At(n+1) = m ;
}
