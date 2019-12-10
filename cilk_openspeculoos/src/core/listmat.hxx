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

#ifndef LISTMAT
#define LISTMAT

/*!
 * \file listmat.hxx
 * \brief Manages list of matrices associated to parent elements
 * \author {Lots of authors}
 * \version 1.0
 */

 
class ParentElement ;
class Matrix ;
template <class T> class Vector ;

/*! \class ListMatrices
   * \brief Manages list of matrices associated to parent elements.
   * 
   * A list of matrices is intended to be owned by a parent element. It
   * contains a list ('parents') of parent elements and a list ('matrices')
   * of matrices associated with every parent element.
   * 
   * Every i-th matrix is typically the interpolation matrix (or the 
   * derived-interpolation matrix) from the owner parent to the i-th 
   * parent. 
   * 
   * One of the parents is typically the owner itself; the associated
   * matrix is the interpolation matrix of the owner.
   * 'parents' and 'matrices' are of the same size; they grow together as
   * additional parent/matrix pairs are accommodated.
   * 
   * FIXME:IMPROVEMENT: ALLOWING DIAGONAL MATRICES TO BE STORED!!!
   */ 
class ListMatrices
{  

protected :
  
    Vector<ParentElement*>* parents ;
    Vector<Matrix*>*        matrices ;

public :

    ListMatrices () ;
    ~ListMatrices () ;

    Matrix* Search (ParentElement* p) ;
    void    Store  (ParentElement* p, Matrix* m) ;
} ;


#endif
