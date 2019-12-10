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

#ifndef VERTEX
#define VERTEX

#include "core/element.hxx"

/*!
 * \file vertex.hxx
 * \brief Manages spectral elements of dimension 0
 * \author {Lots of authors}
 * \version 1.0
 */

/*! \class Vertex
   * \brief Manages spectral elements of dimension 0.
   * 
   * The position of a vertex is characterized by its <{point}>.
   * 
   * A vertex has two main purposes:
   * - to define 1D elements (edges);
   * - to hold values, like any element. Vertices, as parts of 0D meshes,
   *   are heavily used in mortar fields for ensuring C0 inter-element
   *   continuity (with or without mortars).
   * 
   * The attributes 'distributions' and 'ownsDistributions' inherited from
   * Element are not used.
   * 
   * Vertices are usually created by the user (using points) and by the
   * edges, faces and volumes (in the generation methods).
   */ 

class Vertex : public Element
{

protected :

    Point* point ;

public :

    // = CONSTRUCTORS

    Vertex (Point* p) ;
    ~Vertex () ;

    // = TOPOLOGICAL OPERATIONS

    Point*        GetBarycenter () ;
    inline int    GetDimension () ;
    inline Point* GetPoint () ;
    boolean       IsIdenticalTo (Element* elem) ;
    boolean       IsVertex () ;
    void          Print () ;
    int           SearchLocationOf (Element* elem) ;

    inline Vertex*  GetVertex (int i){ return this;}
    inline int  GetNbVertices ( ){ return IONE;}
    // = COORDINATE FIELD

    void          SetToCoordinates (int index) ;

    // = COPY, ADD SPECTRAL ELEMENT VALUES

    void          CopyAddValuesFromVertex (Vertex* elemX, int indX, int indY,
                                           int icx, int icy, CopyAddMode cam,
                                           InterpMode im, NonConform nc,
                                           ReceiveMode rm) ;
    void          CopyAddValuesFromEdge   (Edge* elemX, int indX, int indY,
                                           int icx, int icy, CopyAddMode cam,
                                           InterpMode im, NonConform nc,
                                           ReceiveMode rm) ;
    void          CopyAddValuesFromFace   (Face* elemX, int indX, int indY,
                                           int icx, int icy, CopyAddMode cam,
                                           InterpMode im, NonConform nc,
                                           ReceiveMode rm) ;
} ;


//__________________________ Vertex (inline) __________________________________


int Vertex :: GetDimension () 
{
    // Returns 0.
 
    return 0 ; 
}


Point* Vertex :: GetPoint () 
{
    // Returns the geometrical point of the receiver.
 
    return point ; 
}


#endif
