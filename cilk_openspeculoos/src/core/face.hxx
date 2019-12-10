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

#ifndef FACE
#define FACE

/*!
 * \file face.hxx
 * \brief Management of the faces
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/volume.hxx"
#include "core/edge.hxx"
#include "core/vertex.hxx"
#include "core/distrib.hxx"
#include "core/mesh.hxx"
#include "core/point.hxx"
#include "core/parent.hxx"
#include "core/operator.hxx"
#include "core/matrix.hxx"
#include "core/util.hxx"
#include "core/param.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <math.h>
#include <stdio.h>
#include "core/funct.hxx"
#include "core/element.hxx"

class Mesh2D ;
class OperatorQ1D ;

enum { SOUTH=1, NORTH, WEST, EAST } ;


/*! \class Face
   * \brief Spectral element of dimension 2.
   * 
   * A face is characterized by a set of edges, typically 2 or 3,  stored
   * in <{edges}>. An orientation (positive or negative) is usually 
   * associated by the face to each of its edges (attribute
   * <{orientations}>).
   * 
   * Two adjacent edges of a face point to the same vertex, not just to two
   * geometrically coincident vertices. Similarly, two adjacent faces in a
   * mesh point to a same edge, not just to two coincident edges.
   * 
   * The edges of a face may be straight or curved; the overall surface of
   * a face can be planar or warped. Faces can be used in 2D or 3D 
   * problems.
   * 
   * Currently there is only one subclass of Face: Quad. Maybe these two
   * classes can be shrunk into one single class. 
   */ 

class Face : public Element
{ 

protected :

    Vector<Edge*>* edges ;
    Vector<int>*   orientations ;

    Face () ;

public :

    // = CONSTRUCTORS

    virtual ~Face () ;

    // = MANAGEMENT OF THE EDGES

    virtual int     GetDirectionOf (Edge* e) = 0 ;
    Edge*           GetEdge (int i) ;
    inline int      GetNbSubelements () ;
    virtual int     GetNbEdges () = 0 ;
    Element*        GetSubelement (int i) ;
    boolean         Has (Element* elem) ;
    void            Substitute (Edge* source, Edge* target) ;
    void            SubstituteCommonEdges (Face* face) ;

    // OTHER TOPOLOGICAL OPERATIONS

    virtual Mesh2D* GenerateMesh () = 0 ;
    inline int      GetDimension () ; 
    int             SearchLocationOf (Element* edge) ;

    // = COMPUTATION OF A POINT

    virtual Point*  Eval (real r, real s) = 0 ;

    // = COMPARISON

    virtual boolean IsIdenticalTo (Element* elem) = 0 ;
    boolean         IsFace () ;
    virtual boolean IsQuad () ;

    // = OPERATIONS ON VALUES

    void            GetSetInteriorDof (int index, int iop, int idof, 
                                       real& val) = 0 ;
    // = MORTAR OPERATIONS

    virtual void    ApplyQ (ElementaryField* border, ElementaryField* mortar, 
                            int icb, int icm, InterpMode im) = 0 ;

    // = PRINTING AND DEBUGGING

    virtual void    Print () = 0 ;
    void            PrintPointers () ;
    Point*          GetBarycenter () ;
} ;


/*! \class Quad
   * \brief Quadrilateral element of dimension 2.
   * 
   * A quad is characterized by a set of 4 edges and their orientation.
   * By convention, edge 1 is referred to as the South edge,
   *                edge 2 is referred to as the North edge,
   *                edge 3 is referred to as the West  edge,
   *                edge 4 is referred to as the East  edge.
   * The {S,N,W,E} nomenclature was preferred to {bottom,top,left,right}
   * nomenclature in order avoid confusions with 3D elements (hexahedrons).
   * 
   * The edges are oriented:
   * - the South edge is said to have a positive orientation (i.e., 
   *   orientations(1)=+1) if its first vertex also belongs to the West
   *   edge (rather than to the East edge);
   * - same convention for the North edge;
   * - the West edge is oriented positively if its first vertex also 
   *   belongs to the South edge;
   * - same convention for the East edge.
   *  
   * To summarize, the positive convention is:
   *                South and North edges: from West  to East;
   *                West  and East  edges: from South to North.
   * Note that the edges may or may not be straight.
   * 
   * Any point P(x,y,z) in the face can be obtained by a mapping F(r,s) 
   * with r and s in [-1,1]. This function F is obtained by interpolation
   * of its values on the 4 edges. The Gordon-Hall linear transfinite-
   * interpolation method is used:
   * F(r,s) = phi1(r) F(-1, s) + phi2(r) F(1,s) + phi1(s) F( r,-1) + phi2(s) F(r,1) - phi1(r) phi1(s) F(-1,-1) - phi1(r) phi2(s) F(-1, 1) - phi2(r) phi1(s) F( 1,-1) - phi2(r) phi2(s) F( 1, 1)
   * 
   * where:
   * - the parameters r and s are the supports of the curvilinear 
   *   abscissae in the 2 directions; they both range in [-1,1];
   *   r=-1 -> West  edge;    r=1 -> East  edge;
   *   s=-1 -> South edge;    s=1 -> North edge;
   * - phi1(r) = (1-r)/2 and phi2(r) = (1+r)/2;
   * - F(-1,s), F(1,s), F(r,-1) and F(r,1) are the parametric 
   *   representations of the 4 edges: West, East, South, North, 
   *   respectively;
   * - F(-1,-1), ..., F(1,1) are the 4 corner points.
   * With respect to Gordon-Hall's original formulation, a minor 
   * modification is introduced: the curvilinear abscissae r and s range in
   * [-1,1], rather than in [0,1], in order to gain consistency with class
   * Edge, which is defined on [-1,1].
   * Reference: William J. Gordon and Charles A. Hall, Construction of 
   * Curvilinear Co-ordinate Systems and Application to Mesh 
   * Generation, International Journal for Numerical Methods in 
   * Engineering, Vol 7, pp 461-477, 1973.
   * 
   * A quad is intended to be used either as the geometrical support of
   * one (or possibly several) spectral elements, or as a mesh generator.
   * In the latter case, attribute <{distributions}> contains 2 
   * distributions, one for the West-East direction, one for the South-
   * North direction.
   * 
   * If the quad is used as a generation element, and if one (or more) of
   * its edges is curved (i.e., is not a Line), then the edges generated
   * inside the quad are:
   * - curved (i.e., of type TransfiniteEdge) if the class attribute
   *   'curvedInteriorEdges' is set to 'true'; 
   * - straight (i.e., of type Line) if 'curvedInteriorEdges' is set to 
   * 'false'.
   * The former is the default choice, because in a 3D case if the face is
   * a curved face located on the mesh boundary, then it is important that
   * the geometry of the defined sub-edges of that face be modelled very
   * precisely (example: flow around a cylinder).
   * Otherwise, and typically in 2D cases, the latter is advisable,
   * because computations on straight geometries are more efficient and
   * more accurate.
   * Attribute 'curvedInteriorEdges' can be modified at any time by means
   * of the class method 'SetCurvedInteriorEdges'.
   * 
   * A quad is typically created either by the user (by means of a set of
   * 4 edges) or by a quad used as a mesh generator, see method 
   * 'GenerateMesh'.
   */ 

class Quad : public Face
{ 

protected :

    static int     numberOfEdges ;
    static boolean curvedInteriorEdges ;

    boolean        HasConsistentEdges () ;

public :

    // = CONSTRUCTORS

    Quad (Edge* s, Edge* n, Edge* w, Edge* e) ;
    ~Quad () ;

    // = TOPOLOGICAL OPERATIONS

    int         GetDirectionOf (Edge* e) ;
    inline int  GetNbEdges () ;
    inline int  GetNbVertices () ;
    int         GetOrientationOfEdge (int i) ;
    real        GetTypicalLength (int dir) ;
    Vertex*     GetVertex (int i) ;
    boolean     IsInItsSpace () ;
    void        Print () ;

    // = COORDINATE FIELD

    void        SetToCoordinates (int index) ;

    // EVALUATION OF F(r,s)

    Point*      Eval (real r, real s) ;

    // GENERATION

    Mesh2D*     GenerateMesh () ;
    Edge*       ExtractEdgeBetween (real r1, real r2, real s1, real s2) ;
    static void SetCurvedInteriorEdges (boolean choice) ;

    // = COMPARISON

    boolean     IsIdenticalTo (Element* elem) ;
    boolean     IsQuad () ;

    // = COPY, ADD SPECTRAL ELEMENT VALUES

    void        CopyAddValuesFromEdge   (Edge* elemX, int indX, int indY,
                                         int icx, int icy, CopyAddMode cam,
                                         InterpMode im, NonConform nc,
                                         ReceiveMode rm) ;
    void        CopyAddValuesFromFace   (Face* elemX, int indX, int indY,
                                         int icx, int icy, CopyAddMode cam,
                                         InterpMode im, NonConform nc,
                                         ReceiveMode rm) ;
    void        CopyAddValuesFromVolume (Volume* elemX, int indX, int indY,
                                         int icx, int icy, CopyAddMode cam,
                                         InterpMode im, NonConform nc,
                                         ReceiveMode rm) ;

    // = MORTAR OPERATIONS

    void        ApplyQ (ElementaryField* border, ElementaryField* mortar, 
                        int icb, int icm, InterpMode im) ;

    // = OPERATIONS ON VALUES

    void        GetSetInteriorDof (int index, int iop, int idof, real& val) ;
  
} ;


// ______________________________ Face (inline) _______________________________


int Face :: GetDimension () 
{
    // Returns 2.
 
    return 2 ; 
} 


int Face :: GetNbSubelements ()
{
    // Returns the number of edges of the receiver.

    return GetNbEdges() ;
}


//_______________________________ Quad (inline) _______________________________


int Quad :: GetNbEdges ()
{
    // Returns the number of edges of the receiver: 4.

//  return numberOfEdges ;        // issues core dump
    return 4 ;
}


int Quad :: GetNbVertices ()
{
    // Returns the number of vertices of the receiver: 4.

    return 4 ;
}


#endif
