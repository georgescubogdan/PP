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

#ifndef EDGE
#define EDGE

/*!
 * \file edge.hxx
 * \brief Spectral element
 * \author {Lots of authors}
 * \version 1.0
 */


using namespace std ;

#include "core/face.hxx"
#include "core/mesh.hxx"
#include "core/point.hxx"
#include "core/parent.hxx"
#include "core/distrib.hxx"
#include "core/operator.hxx"
#include "core/matrix.hxx"
#include "core/util.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <math.h>
#include <stdio.h>
#include "core/element.hxx"
#include "core/vertex.hxx"


class Mesh1D ;
class Quad ;
class Vertex ;
class OperatorQ1D ;


/*! \class Edge
   * \brief Spectral element of dimension 1.
   * 
   * An edge is characterized by two vertices 'vertex1' and 'vertex2'.
   * The edge's shape between these vertices is defined in the derived
   * classes below.
   * When two edges have a geometrically coincident vertex, this vertex
   * must be the same object, in order to make subsequent topological
   * operations successful.
   * 
   * If an edge is used as a generation object, its attribute 
   * <{distributions}> contains 1 distribution.
   * 
   * Edges are created by means of vertices, typically by the user or by
   * the generation methods of classes Edge, Face and Volume.
   */
class Edge : public Element
{

protected :

    Vertex*               vertex1 ;
    Vertex*               vertex2 ;
    Vector<OperatorQ1D*>* operatorsQ ;

    Edge () ;

public :

    // = CONSTRUCTORS

    Edge (Vertex* v1, Vertex* v2) ;
    virtual ~Edge () ; 

    // = TOPOLOGICAL OPERATIONS

    Mesh1D*         GenerateMesh () ;
    inline int      GetDimension () ;
    boolean         HasSameOrientationAs (Edge* edge) ;

    // = OPERATIONS ON THE 2 VERTICES

    inline int      GetNbSubelements () ;
    inline int      GetNbVertices () ;
    Element*        GetSubelement (int i) ;
    inline Vertex*  GetVertex1 () ;
    inline Vertex*  GetVertex2 () ;
    Vertex*         GetVertex (int i) ;
    Vertex*         GetVertexAt (Point* p) ;
    Vertex*         GetOtherVertex (Vertex* v) ;
    inline boolean  Has (Element* elem) ;
    boolean         HasVertexAt (Point* p) ;
    int             SearchLocationOf (Element* v) ;
    void            SetVertex (int i, Vertex* v) ;
    void            SetVerticesTo (Vertex* v1, Vertex* v2) ;
    void            Substitute (Vertex* source, Vertex* target) ;
    void            SubstituteCommonVertices (Edge* edge) ;
 
    // COMPUTATION OF A POINT

    virtual Point*  Eval (real r) = 0 ;

    // = COMPARISON, DUPLICATION

    virtual Edge*   DuplicateBetween (real r1, real r2) = 0 ;
    virtual boolean IsIdenticalTo (Element* elem) = 0 ;
    inline  boolean IsEdge () ;
    virtual boolean IsCircle () ;
    virtual boolean IsLine () ;
    virtual boolean IsTransfiniteEdge () ;

    // = COORDINATE FIELD

    void            SetToCoordinates (int index)  ;
 
    // = COPY, ADD SPECTRAL ELEMENT VALUES

    void            CopyAddValuesFromVertex (Vertex* elemX, int indX, int indY,
					     int icx, int icy, CopyAddMode cam,
					     InterpMode im, NonConform nc,
					     ReceiveMode rm) ;
    void            CopyAddValuesFromEdge (Edge* elemX, int indX, int indY,
					   int icx, int icy, CopyAddMode cam,
					   InterpMode im, NonConform nc,
					   ReceiveMode rm) ;
    void            CopyAddValuesFromFace (Face* elemX, int indX, int indY,
					   int icx, int icy, CopyAddMode cam,
					   InterpMode im, NonConform nc,
					   ReceiveMode rm) ;
    // = OPERATIONS ON VALUES

    void            GetSetInteriorDof (int index, int iop, int idof, 
                                       real& val) ;

    // = MORTAR CONSTRAINTS

    FullMatrix*     GetMatrixQ (ParentElement* border, ParentElement* mortar) ;

    // = MISCELLANEOUS

    Point*          GetBarycenter () ;
    virtual void    Print () ;
    void            PrintPointers () ;
} ;


/*! \class Line
   * \brief Straight geometrical element of dimension 1.
   * 
   * [Vincent Keller] FIXME : No documentation 
   */

class Line : public Edge
{

public :

    Line (Vertex*, Vertex*) ;
    ~Line () ;

    Edge*   DuplicateBetween (real r1, real r2) ;
    Point*  Eval (real r) ;
    real    GetLength () ;
    real    GetTypicalLength (int dir) ;
    boolean IsIdenticalTo (Element* elem) ;
    boolean IsLine () ;
    void    Print () ;
} ;



/*! \class Circle
   * \brief Geometrical element of dimension 1 with arc-circular shape.
   * 
   * A circular edge in a 2D or 3D space is defined by a <{center}> and a
   * <{radius}>, besides its two inherited end vertices.
   * 3D case is limited to the case where the two vertices and the center
   * lie in the same horizontal plane Z=cst.
   * If the arc encompasses a complete circle, the two end vertices are
   * identical.
   * Like other curves, the curvilinear abscissa ranges in [-1,1] ; in
   * order to speed up the calculations of points (method "Eval"), the
   * radian angles <{angle1}> and <{angle2}> (in the plane defined by the
   * arc and the center) of the two vertices are stored.
   */ 

class Circle : public Edge
{

protected :

    Point* center ;
    real   radius ;
    real   angle1 ;
    real   angle2 ;

public :

    Circle (Vertex* v1, Vertex* v2, Point* c) ;
    ~Circle () ;

    Edge*         DuplicateBetween (real r1, real r2) ;
    Point*        Eval (real r) ;
    inline Point* GetCenter () ;
    boolean       IsCircle () ;
    boolean       IsIdenticalTo (Element* elem) ;
    void          Print () ;
} ;

/*! \class TransfiniteEdge
   * \brief Geometrical element of dimension 1.
   * 
   * A transfinite edge is defined as the restriction of a general face
   * F(r,s), r and s defined on I=[-1,1] to either
   *   . F(r0,s), r0 given and s ranging on a subinterval [a1,a2] of I, or
   *   . F(r,s0), s0 given and r ranging on a subinterval [a1,a2] of I.
   *    Like edges of other types, from the viewpoint of any client, the
   * curvilinear abscissa (s or r) is defined on [-1,1], not on [a1,a2].
   * 
   * A transfinite edge is characterized by:
   * - <{face}>: a face which defines the F(r,s) function;
   * - <{blockedDirection}>: =1 if r0 is imposed, =2 if s0 is imposed;
   * - <{rs0}>: r0 or s0, depending on <{blockedDirection}>;
   * - {<a1>} and {<a2>}: define the interval [a1,a2]. 
   * 
   * Usually r0 (or s0) is different from -1 and 1, because in those cases
   * the transfinite edge is simply a straight, or circular, etc, edge.
   *     
   * Transfinite edges are typically created in the mesh generation of a
   * general face (not a Rectangle).
   * 
   * Remark: If 2D problems only were relevant, one could dispense with
   * transfinite edges, because the edges inside a face can be appropria-
   * tely generated as straight edges. This does not hold in 3D problems,
   * where faces must be accurately refined, as in the case of generating
   * faces on a cylinder-like shape.
   * Note that this class has no 3D counterpart ("TransfiniteQuad"), becau-
   * se the faces inside a volume can be appropriately generated with
   * straight edges, rather than as restrictions of the volume.
   */

class TransfiniteEdge : public Edge
{

protected :

    Quad* quad ;
    int   blockedDirection ;
    real  rs0 ;
    real  a1 ;
    real  a2 ;

public :

    TransfiniteEdge (Quad* q, int dir, real RS0, real A1, real A2) ;

    Edge*   DuplicateBetween (real rs1, real rs2) ;
    Point*  Eval (real rs) ;
    boolean IsIdenticalTo (Element* elem) ;
    boolean IsTransfiniteEdge () ;
    void    Print () ;
} ;



//______________________________ Edge (inline) ________________________________


int Edge :: GetDimension () 
{
    // Returns 1.

    return 1 ; 
}


int Edge :: GetNbSubelements ()
{
    // Returns the number of vertices of the receiver.

    return 2 ;
}

int Edge :: GetNbVertices ()
{
    // Returns the number of vertices of the receiver.

    return 2 ;
}


Vertex* Edge :: GetVertex1 () 
{
    // Returns the first vertex of the receiver.
 
    return vertex1 ; 
} 


Vertex* Edge :: GetVertex2 ()
{
    // Returns the second vertex of the receiver.

    return vertex2 ; 
}


boolean Edge :: Has (Element* elem)
{
    // Returns True if 'elem' is one of the two vertices of the receiver (same
    // object), else returns False.

    return (elem == vertex1 || elem == vertex2) ;
}


boolean Edge :: IsEdge ()
{
    // Returns True, because the receiver is an edge.

    return true ;
}


//______________________________ Circle (inline) ______________________________


Point* Circle :: GetCenter ()
{
    // Returns the center of the receiver.

    return center ;
}


#endif
