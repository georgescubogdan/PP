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

#ifndef MESH
#define MESH

/*!
 * \file mesh.hxx
 * \brief Manages mesh
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/vector.hxx"
#include "core/garbaged.hxx"
#include <string>
#include <vector>

class Point ;
class Element ;
class Face ;
class Edge ;
class FlatField ;
class ParentElement ;
class Communicator ;


/*! \class Mesh
   * \brief Base class managing meshes of spectral elements of any dimension.
   * 
   * A mesh is primarily a set of spectral elements. All of the spectral
   * elements must be of the same dimension, i.e., their geometry must be
   * described by the same number of local parameters (r,s,t).
   * 
   * A mesh has a <{skeleton}> which contains smaller-dimensional elements.
   * For exemple, a mesh M of 2D elements has a skeleton S of 1D elements.
   * S is made of the edges of all faces of M. The skeleton of S is in turn
   * a mesh of 0D elements. 
   * A mesh is also able to extract its skin (border), which is useful to 
   * build meshes suitable for boundary conditions.
   * 
   * A mesh possesses a set of four "basic" fields, namely: 
   * 
   * - <{coordinates}>: the ndim components of this field correspond to the
   * directions of the mesh, eg, 3 components (X,Y,Z) for a 3D mesh;
   * 
   * - <{jacobian}>: this field has 1 component, the determinant of the
   * partial derivatives of the ndim coordinates (X,Y,Z) with respect to
   * the ndim local coordinates (r,s,t) of every spectral element;
   * 
   * - <{invJacobian}>: the inverse of the jacobian matrix. This field has
   * ndim*ndim components, the partial derivatives of the local 
   * coordinates (r,s,t) with respect to the cartesian coordinates
   * (X,Y,Z);
   * 
   * - <{normal}>: this field has nbCoordinates=ndim+1 components, the
   * components of the normal vector to the mesh.
   * 
   * Remark. A mesh usually as either its <{invJacobian}> or its <{normal}>
   * initialized, but not both:
   * . <{invJacobian}> can be requested only from a "main" mesh, i.e., not
   *   from a mesh created as the skeleton or the skin of another mesh.
   *   Otherwise it has no sense, because the jacobian matrix is not well
   *   defined. For example, the skeleton of a 2D mesh is a 1D mesh (a mesh
   *   made of edges). It has 2 coordinates (X and Y) but every element has
   *   only 1 parameter (r).
   * . <{normal}> can be requested only from a skin mesh, because it has a
   *   meaning only at the boundary.
   * 
   * When the mesh creates its coordinate field, it uses:
   * - <{nbCoordinates}> as the number of components (X,Y,Z). By default,
   * this attribute is set to the dimension of the mesh. This value may
   * need to be overridden. For example, the skin (or the skeleton) of a
   * 2D mesh is a 1D mesh. For that mesh the default value (1) of
   * <{nbCoordinates}> must be overridden and set to 2.
   * 
   * - for every spectral element, in every direction, the interpolation 
   * scheme given in attribute <{defaultParentElement}>. By default, this
   * attribute is set to GLL of degree 1. The user can override this
   * default value.
   * 
   * When the mesh creates its <{jacobian}> determinant, its
   * <{invJacobian}> matrix or its <{normal}> vector, it uses the same 
   * interpolation as for its <{coordinates}> field.
   * 
   * Note that whatever the type and the degree of the chosen default
   * parent element may be, the values of the 4 basic fields are exact,
   * because the definition of every element's geometry is analytic.
   * 
   * A mesh is able to "dispatch" its elements, i.e., to assign each of 
   * them to a processor. This assignment operation:
   * - must be done before a field is constructed on the mesh, so that the
   * elements know where to store the values of their elementary fields;
   * - need not be done on a scalar implementation, because in that case
   * all elements are assigned to the same processor 0;
   * - can be customized by the user, possibly resulting in a better load
   * balancing and/or a smaller amount of processor communication
   * (e.g., the general-purpose method DispatchElements vs. the 
   * customized method DispatchElementsByBlocks).
   */ 

class Mesh : public Vector<Element*>, public GarbagedObject
{

protected :

    Mesh*                  skeleton ;                    // skeleton mesh

    int                    nbCoordinates ;               // basic fields
    FlatField*             coordinates ;
    FlatField*             jacobian ;
    FlatField*             invJacobian;
    FlatField*             normal ;
    ParentElement*         defaultParentElement ;
    boolean                inhibitedCoordinates ;

    Vector<Communicator*>* internalCommunicators ;       // communicators
    Vector<Communicator*>* sendReceiveCommunicators ;
    Vector<Communicator*>* receiveCommunicators ;

public :

    // = CONSTRUCTORS

    Mesh () ;
    virtual ~Mesh () ;

    // = TOPOLOGICAL OPERATIONS

    static Mesh*   CreateMesh (int dim) ;
    Mesh*          Except (Mesh* m) ;
    Mesh*          Extract (int dir, real c) ;
    virtual int    GetDimension () ;
    Element*       GetElementCloserTo (Point* p) ;
    Mesh*          GetElementsEncompassing (Element* elem) ;
    virtual Mesh*  GetSkeleton () ;
    virtual Mesh*  GetSkin () ;
    virtual void   Put (Element* elem) = 0 ;

    // = MANAGEMENT OF THE BASIC FIELDS (COORDINATES, JACOBIAN, INVERSE
    //   JACOBIAN, NORMAL)

    void           DeleteBasicFields () ;
    void           DeleteJacobian () ;
    void           DeleteInvJacobian () ;
    FlatField*     GetCoordinates () ;
    FlatField*     GetInvJacobian () ;
    FlatField*     GetJacobian () ;
    int            GetNbCoordinates () ;
    FlatField*     GetNormal (RealVector* external) ;
    void           InhibitCoordinates () ;
    void           InterpolateCoordinatesLike (Mesh* m) ;
    void           SetDefaultParentElement (ParentElement* p) ;

    // = PARALLEL PROCESSING

    void           DispatchElements () ;
    Communicator*  GetInternalCommunicatorTo (Mesh* m) ;
    Communicator*  GetReceiveCommunicatorTo (Mesh* m) ;
    Communicator*  GetSendReceiveCommunicatorTo (Mesh* m) ;
    boolean        IsDispatched () ;
    void           PrintDispatching () ;
    // = DEBUGGING

    void           Print () ;
    void           PrintPointers () ;
    void           PrintBarycenters () ;
    void           PrintConnectivity () ;


    // = GARBAGE COLLECTION

    void           BreakCyclicReferences () ;
    std::vector<int> const& GetLocalElements() const;
protected:
    void ComputeLocalElements();
private:
    std::vector<int> localElements;
} ;


/*! \class Mesh0D
   * \brief Manages meshes made of vertices.
   * 
   * Manages meshes made of vertices.
   */ 
class Mesh0D : public Mesh
{

public :
  
    // = CONSTRUCTORS
  
    Mesh0D () ;
    ~Mesh0D () ;
  
    // = TOPOLOGICAL OPERATIONS
  
    void       Put (Element* elem) ;
    inline int GetDimension () ;
    Mesh*      GetSkeleton () ;
    Mesh*      GetSkin () ;
  
} ;



/*! \class Mesh1D
   * \brief Manages meshes made of edges.
   * 
   * Manages meshes made of edges.
   */ 
class Mesh1D : public Mesh
{

public :

    // = CONSTRUCTORS

    Mesh1D () ;
    ~Mesh1D () ;

    // = TOPOLOGICAL OPERATIONS

    Mesh1D*     ExtractMeshBetween (Point* pA, Point* pB, Point* pC) ;
    Mesh1D*     ExtractMeshBetween (Point* pA, Point* pB) ;
    inline int  GetDimension () ;
    Mesh1D*     GetEdgesWith (Point* p) ;
    void        Put (Element* elem) ;
    void        Put (Mesh1D* mesh2) ;
} ;



/*! \class Mesh2D
   * \brief Manages meshes made of faces
   * 
   * Manages meshes made of faces.
   */ 
class Mesh2D : public Mesh
{

public :

    // = CONSTRUCTORS

    Mesh2D () ;
    ~Mesh2D () ;
  
    // = TOPOLOGICAL OPERATIONS

    inline int  GetDimension () ;
    void        Put (Element* elem) ;
    void        Put (Mesh2D* mesh2) ;
 
    // = PARALLEL PROCESSING

    void        DispatchElementsByBlocks (int nbElem1, int nbElem2) ;
} ;



/*! \class Mesh3D
   * \brief Manages meshes made of volumes
   * 
   * Manages meshes made of volumes.
   */ 
class Mesh3D : public Mesh
{

public :

    // = CONSTRUCTORS

    Mesh3D () ;

    // = TOPOLOGICAL OPERATIONS

    inline int  GetDimension () ;
    void        Put (Element* elem) ;
    void        PutWithoutCheck (Element* elem) ;
    void        Put (Mesh3D* mesh2) ;
    // = PARALLEL PROCESSING

    void        DispatchElementsByBlocks (int nbElem1, int nbElem2, 
                                          int nbElem3) ;
} ;


//______________________________ Mesh0D (inline) ______________________________


int Mesh0D :: GetDimension () 
{
    // Returns the geometrical dimension of the receiver: 0 .

    return 0 ; 
} 

//______________________________ Mesh1D (inline) ______________________________


int Mesh1D :: GetDimension () 
{
    // Returns the geometrical dimension of the receiver: 1.

    return 1 ; 
} 


//______________________________ Mesh2D (inline) ______________________________


int Mesh2D :: GetDimension () 
{
    // Returns the geometrical dimension of the receiver: 2.

    return 2 ; 
}


//______________________________ Mesh3D (inline) ______________________________


int Mesh3D :: GetDimension () 
{
    // Returns the geometrical dimension of the receiver: 3.

    return 3 ; 
}


#endif
