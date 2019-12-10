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

#ifndef OPERATORQ
#define OPERATORQ

/*!
 * \file operators.hxx
 * \brief Manages operators
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/matrix.hxx"

class ParentElement ;


/*! \class OperatorQ1D
   * \brief FIXME : no doc
   * 
   * FIXME [Vincent Keller] No documentation !
   */ 
class OperatorQ1D
{
protected :
    
    ParentElement* border ;
    ParentElement* mortar ;
    FullMatrix*    matrix;

public :

    inline OperatorQ1D (ParentElement* b, ParentElement* m, FullMatrix* q) ;
    inline ~OperatorQ1D () ;

    inline ParentElement* GetBorderParent () ;
    inline ParentElement* GetMortarParent () ;
    inline FullMatrix*    GetMatrix () ;
} ;


//____________________________ OperatorQ1D (inline) ___________________________


OperatorQ1D :: OperatorQ1D (ParentElement* b, ParentElement* m, FullMatrix* q)
{
    // Constructor.

    border = b ;
    mortar = m ;
    matrix = q ;
}


OperatorQ1D :: ~OperatorQ1D ()
{
    // Destructor.

    delete matrix ;
}


ParentElement* OperatorQ1D :: GetBorderParent ()
{
    return border ;
}


ParentElement* OperatorQ1D :: GetMortarParent ()
{
    return mortar ;
}


FullMatrix* OperatorQ1D :: GetMatrix ()
{
    return matrix ;
}


#endif
