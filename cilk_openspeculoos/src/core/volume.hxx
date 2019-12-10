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

#ifndef VOLUME
#define VOLUME

/*!
 * \file volume.hxx
 * \brief Manages spectral elements of dimension 3
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/element.hxx"
#include "core/vector.hxx"
#include <string>

class Mesh3D ;
class Quad ;

/*! \class Volume
   * \brief Manages specral elements of dimension 3.
   * 
   * Currently there is only one subclass of Volume: Hexa. Maybe these two
   * classes can be shrunk into one single class.
   */ 

class Volume : public Element
{

public :

    // = CONSTRUCTORS

    Volume () ;
    virtual ~Volume () ;

    // = TOPOLOGICAL OPERATIONS
 
    virtual Point*       Eval (real r, real s, real t) = 0 ;
    virtual Mesh3D*      GenerateMesh () = 0 ;
    inline int           GetDimension () ; 
    virtual void         GetDirectionsOf (Face* f, int& dir1, int& dir2) = 0 ;
    virtual void         GetDirectionsOfFace (int i, int& dir1, int& dir2) = 0;
    virtual Face*        GetFace (int i) = 0 ;
    virtual int          GetNbFaces () = 0 ;
    virtual int          GetNbEdges () = 0 ;
    inline int           GetNbSubelements () ;
    Element*             GetSubelement (int i) ;

    virtual boolean      Has (Element* elem) = 0 ;
    inline int           IsVolume () ;
    virtual int          SearchLocationOf (Element* face) = 0 ;
    virtual void         Substitute (Face* source, Face* target) = 0 ;

    // = PRINTING AND DEBUGGING

    Point*               GetBarycenter () ;
    virtual void         Print () = 0 ;
    void                 PrintPointers () ;
 
    // = OPERATIONS ON FIELDS

    virtual void             SetToCoordinates (int index) = 0 ;
    virtual ElementaryField* Extract (int indy, int icy, Face* face,
                                      ReceiveMode rm) = 0 ;

    // = OPERATIONS ON VALUES

    void                 GetSetInteriorDof (int index, int iop, int idof,
                                            real& val) = 0 ;
} ;


/*! \class Hexa
   * \brief 6-face geometrical element of dimension 3.
   * 
   * A hexahedron is characterized by a set of 6 quadrilateral faces and
   * their orientation. 
   * By convention, face 1 is referred to as the left   face,
   *                face 2 is referred to as the right  face,
   *                face 3 is referred to as the front  face,
   *                face 4 is referred to as the rear   face,
   *                face 5 is referred to as the bottom face,
   *                face 6 is referred to as the top    face.
   * 
   * The hexahedron associates an orientation to each face. The orientation
   * of every face is defined as a (rotation,yaw) pair. For every face
   * there are 8 such possible orientation pairs.
   * 
   * The possible values taken by the rotation index of every face are 0, 
   * 1, 2 and 3. The convention for a rotation index equal to 0 is:
   *    left   face: its South edge also belongs to the bottom face;
   *    right  face: its South edge also belongs to the bottom face;
   *    front  face: its South edge also belongs to the bottom face;
   *    rear   face: its South edge also belongs to the bottom face;
   *    bottom face: its South edge also belongs to the front  face;
   *    top    face: its South edge also belongs to the front  face.
   * A rotation index equal to 1, 2 or 3 means that, with respect to a 
   * rotation index equal to 0, the face is rotated 90, 180 or 270 degrees
   * counterclockwise, respectively, looking to the face from OUTSIDE the
   * volume. For example, for the left face this yields:
   *    rot = 0: its South edge also belongs to the bottom face;
   *    rot = 1: its South edge also belongs to the front  face;
   *    rot = 2: its South edge also belongs to the top    face;
   *    rot = 3: its South edge also belongs to the rear   face.
   * 
   * The possible values taken by the yaw index are 1 and 2. The convention
   * for a yaw index equal to 1 is: when the face is looked at from OUTSIDE
   * the volume, its 3rd (West) edge stands on the left and its 4th (East)
   * edge stands on the right.
   * The yaw index is equal to 2 if this definition holds if the face is
   * looked at from INSIDE the volume.
   * Equivalent definition: the yaw index is set to 1 if: let 'r' be the 
   * support of the curvilinear abscissa from edge 3 (West) to edge 4 
   * (East) and 's' be the support of the curvilinear abscissa from edge 1 
   * (South) to edge 2 (North); then the cross product 't' of r and s 
   * points OUTSIDE the volume.
   * 
   * Note that on top of this, every face associates an orientation to its
   * edges, every edge having 2 possible orientations.
   * 
   * The faces of the hexahedron can have any shape. In particular they can
   * be warped (ie, not in a plane). 
   * Any point P(x,y,z) in the volume can be obtained by a mapping F(r,s,t) 
   * with r, s, and t in [-1,1]. This function F is obtained by 
   * interpolation of its values on the 6 faces. The Gordon-Hall linear 
   * transfinite-interpolation method is used:
   * F(r,s,t) = phi1(r) F(-1, s, t) + phi2(r) F(1,s,t)
   *          + phi1(s) F( r,-1, t) + phi2(s) F(r,1,t)
   *          + phi1(t) F( r, s,-1) + phi2(s) F(r,s,1)
   *          - phi1(r) phi1(s) F(-1,-1, t) - phi1(r) phi2(s) F(-1, 1, t)
   *          - phi2(r) phi1(s) F( 1,-1, t) - phi2(r) phi2(s) F( 1, 1, t)
   *          - phi1(s) phi1(t) F( r,-1,-1) - phi1(s) phi2(t) F( r,-1, 1)
   *          - phi2(s) phi1(t) F( r, 1,-1) - phi2(s) phi2(t) F( r, 1, 1)
   *          - phi1(t) phi1(r) F(-1, s,-1) - phi1(t) phi2(r) F( 1, s,-1)
   *          - phi2(t) phi1(r) F(-1, s, 1) - phi2(t) phi2(r) F( 1, s, 1)
   *          + phi1(r) phi1(s) phi1(t) F(-1,-1,-1)
   *          + phi1(r) phi1(s) phi2(t) F(-1,-1, 1)
   *          + phi1(r) phi2(s) phi1(t) F(-1, 1,-1)
   *          + phi1(r) phi2(s) phi2(t) F(-1, 1, 1)
   *          + phi2(r) phi1(s) phi1(t) F( 1,-1,-1)
   *          + phi2(r) phi1(s) phi2(t) F( 1,-1, 1)
   *          + phi2(r) phi2(s) phi1(t) F( 1, 1,-1)
   *          + phi2(r) phi2(s) phi2(t) F( 1, 1, 1)
   * where:
   *   - the parameters r, s and t are the supports of the curvilinear
   *     abscissae in the 3 directions; they all range in [-1,1];
   *       r=-1 -> left face;     r=1 -> right face;
   *       s=-1 -> front face;    s=1 -> rear face;
   *       t=-1 -> bottom face;   t=1 -> top face;
   *   - phi1(x) = (1-x)/2 and phi2(x) = (1+x)/2;
   *   - F(-1,s,t), ..., F(r,s,1) are the parametric representations of the
   *     6 faces: left, right, front, rear, bottom, top, respectively.
   *   - F(-1,-1,t), ..., F(1,s,1) are the parametric representations of
   *     the 12 edges;
   *   - F(-1,-1,-1), ..., F(1,1,1) are the 8 corner points.
   * In other words, r denotes the "left-to-right" abscissa, s the "front-
   * to-rear" abscissa and t the "bottom-to-top" abscissa.
   * Using the directional projectors Pr, Ps and Pt given by Gordon and
   * Hall in the 2D case, the expression for P is:
   *     P = Pr + Ps + Pt - Pr Ps - Ps Pt - Pt Pr + Pr Ps Pt, with
   *     Pr(F) = phi1(r) F(-1,s,t) + phi2(r) F(1,s,t), Ps=..., Pt=... .
   * With respect to Gordon-Hall's original formulation, a minor 
   * modification is introduced: the curvilinear abscissae r, s and t range
   * in [-1,1], rather than in [0,1], in order to gain better compatibility
   * with other classes (Edge, Face) which are defined on [-1,1].
   * Reference: William J. Gordon and Charles A. Hall, Construction of 
   *     Curvilinear Co-ordinate Systems and Application to Mesh 
   *     Generation, International Journal for Numerical Methods in 
   *     Engineering, Vol 7, pp 461-477, 1973.
   * 
   * A hexahedron is typically created either by the user (by means of a
   * set of 6 faces) or by another hexahedron (used as a mesh generator,
   * see method 'GenerateMesh'). When the user creates a hexahedron, he/she
   * does not have to care about rotation/yaw indices of faces: all 
   * combinations are taken into account (and checked). The only 
   * requirement is that the 6 given faces must make a hexahedron; this 
   * implies that, for instance, the faces declared as "front" and "left"
   * have one common edge. Note that these faces must point to the same 
   * edge, rather than just to two geometrically coincident edges.
   * 
   * A hexahedron can be used either as the geometrical support of a 
   * spectral element or as mesh generator. In the former case, the 
   * inherited attribute <{distributions>} is not used; in the latter case 
   * it contains 3 distributions corresponding to the r, s and t directions
   * respectively.
   */ 
class Hexa : public Volume
{ 

protected :

    Vector<Quad*>* faces ;
    Vector<int>*   rotations ;
    Vector<int>*   yaws ;
    static int     numberOfFaces ;

    void CheckConsistentFaces () ;
    void Inject (Face* elemX, ElementaryField* eX, ElementaryField* eY, 
                 int icx, int icy, CopyAddMode cam, ReceiveMode rm) ;
    void ReindexToFace (Face* elemX, ElementaryField* eX, int icx);
    ElementaryField* ReindexToVolume (Face* elemX, ElementaryField* eX, 
				      int icx);
    void VolumeToFace (real r, real s, real t, int f, real& a, real& b) ;
    void VolumeFaceToFace (real A, real B, int f, real& a, real& b);
    void Wrong (int i, string face) ;

public :

    // = CONSTRUCTORS

    Hexa (Quad* le, Quad* ri, Quad* fr, Quad* re, Quad* bo, Quad* to) ;
    ~Hexa () ;

    // = TOPOLOGICAL OPERATIONS

    Mesh3D*          GenerateMesh () ;
    void             GetDirectionsOf (Face* f, int& dir1, int& dir2) ;
    void             GetDirectionsOfFace (int i, int& dir1, int& dir2) ;
    Face*            GetFace (int i) ;
    inline int       GetNbFaces () ;
    inline int       GetNbVertices () ;
    inline int       GetNbEdges () ;
    Vertex*          GetVertex (int i)  ;
    real             GetTypicalLength (int dir) ;
    boolean          Has (Element* elem) ;
    int              SearchLocationOf (Element* face) ;
    void             Substitute (Face* source, Face* target) ;
    enum             { LEFT=1, RIGHT, FRONT, REAR, BOTTOM, TOP } ;

    // POINT EVALUATION

    Point*           Eval (real r, real s, real t) ;

    // = OPERATIONS ON FIELDS

    void             SetToCoordinates (int index) ;
    ElementaryField* Extract (int indy, int icy, Face* face, ReceiveMode rm) ;
    void             CopyAddValuesFromVolume (Volume* elemX, 
					      int indx, int indy, int icx, int icy,
					      CopyAddMode cam, InterpMode im, 
					      NonConform nc, ReceiveMode rm) ;
    void             CopyAddValuesFromFace (Face* elemX, int indX, int indY,
                                            int icx, int icy, CopyAddMode cam,
                                            InterpMode im, NonConform nc,
                                            ReceiveMode rm) ;

    // = OPERATIONS ON VALUES

    void             GetSetInteriorDof (int index, int iop, int idof,
                                        real& val) ;
    // = DEBUGGING

    void             Print () ;
} ;


// _____________________________ Volume (inline) ______________________________


int Volume :: GetDimension () 
{
    // Returns the number of curvilinear abscissae of the receiver: 3.
 
    return 3 ; 
} 


int Volume :: GetNbSubelements ()
{
    // Returns the number of faces of the receiver.

    return GetNbFaces() ;
}


boolean Volume :: IsVolume ()
{
    // Returns True because the receiver is a volume.

    return true ;
}


// ______________________________ Hexa (inline) _______________________________


int Hexa :: GetNbFaces ()
{
    // Returns the number of faces of the receiver: 4.

    return numberOfFaces ;
}
int Hexa :: GetNbVertices ()
{
    // Returns the number of vertices of the receiver

    return 8 ;
}

int Hexa :: GetNbEdges ()
{
    // Returns the number of vertices of the receiver

    return 12 ;
}


#endif
