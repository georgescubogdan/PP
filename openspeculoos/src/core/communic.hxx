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

#ifndef COMMUNIC
#define COMMUNIC

/*!
 * \file communic.hxx
 * \brief Communication within OpenSPECULOOS
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/mesh.hxx"
#include "core/element.hxx"
#include "core/point.hxx"
#include "core/util.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include "core/vector.hxx"

#include <stdio.h>
#include <stdlib.h>

class Mesh ;
class Element ;

enum {INTERNAL_COMMUNICATOR, SEND_RECEIVE_COMMUNICATOR, RECEIVE_COMMUNICATOR} ;

/*! \class Communicator
   * \brief Manages a list of pairs of elements assigned to different processors.
   *
   * 
   * 
   * A communicator is intended to improve efficiency on a parallel 
   * implementation when a field defined on a mesh 1 are copied on (or 
   * added to) a field defined on another mesh 2. Its effect is
   * particularly noticeable in the assembly and distribution operations.
   * 
   * A communicator is characterized by two meshes <{mesh1}> and <{mesh2}>.
   * It consists primarily in <{pairs}> which is a list of pairs 
   * elem1/elem2, where elem1 belongs to <{mesh1}> and elem2 to <{mesh2}>.
   * 
   * When a field on mesh 1 must send its values to a field on mesh 2
   * (e.g., in an assembly operation mesh 1 can be a 2D mesh and mesh 2 its
   * 1D skeleton), on every processor the field directly gets the pairs
   * elem1/elem2 relevant to that processor, instead of scanning all 
   * elements of the meshes to search for relevant pairs. Avoiding such
   * long search loops strongly improves computational efficiency.
   *     
   * Communicators come in 3 flavours, denoted by <{type}>: 
   * internal communicator: in any pair elem1/elem2, both elem1 and elem2
   * are assigned to the current processor. No communication will take
   * place;
   * send-receive communicator: in any pair, either elem1 or elem2 (but
   * not both) is assigned to the current processor. If, say, elem1 is
   * assigned to the current processor and elem2 to another processor p, 
   * then:
   * - the current processor will send values to the other processor
   * during the next communication phase;
   * - processor p will know that, during the next communication phase,
   * it will receive values from the current processor. It will also
   * how many values it will receive;
   * receive communicator: elem2 (but not elem1) is on the current
   * processor. During the next communication phase, the current 
   * processor will receive the data from the processor of elem1.
   * 
   * On the reason why the 2nd type above is "send-receive", and not just
   * "send", see file parallel.hxx.
   * 
   * In a serial implementation, only communicators of type "internal" are
   * relevant. For example, if <{mesh1}> is a 2D mesh with n elements and
   * <{mesh2}> is its 1D skeleton, the communicator will consist of the 4*n
   * pairs elem1/elem2 where elem2 is an edge of elem1.
   * 
   * Remarks:
   * - On two different processors, the communicators are different. This
   * is in contrast with the meshes, which are identical on all 
   * processors: all processors store all elements;
   * - <{mesh1}> and <{mesh2}> can be of the same dimension. They can even
   * be the same mesh, in which case any pair turns out to contain twice
   * the same element.
   */

class Communicator
{

protected :

    int                type ;
    Mesh*              mesh1 ;
    Mesh*              mesh2 ;
    Vector<Element**>* pairs ;

    void CreatePairs () ;

public :

    // = CONSTRUCTORS

    Communicator (Mesh* m1, Mesh* m2, int type) ;
    ~Communicator () ;

    // = INQUIRY

    inline Mesh*     GetMesh1 () ;
    inline Mesh*     GetMesh2 () ;
    inline Element** GetPair (int i) ;
    inline int       GetSize() ;

    // = DEBUGGING

    void             Print (int p=-1) ;
    void             PrintPair (int i) ;
} ;


//___________________________ Communicator (inline) ___________________________


Element** Communicator :: GetPair (int i)
{
    // Returns the 'i'-th pair elem1/elem2 of the receiver.

    return pairs->At(i) ;
}


Mesh* Communicator :: GetMesh1 ()
{
    // Returns the first of the two meshes involved in the receiver.

    return mesh1 ;
}


Mesh* Communicator :: GetMesh2 ()
{
    // Returns the second of the two meshes involved in the receiver.

    return mesh2 ;
}


int Communicator :: GetSize ()
{
    // Returns the number of pairs elem1/elem2 of the receiver.

    return pairs->GetSize() ;
}


#endif
