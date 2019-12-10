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

// mesh.cxx

#include "core/mesh.hxx"
#include "core/point.hxx"
#include "core/vertex.hxx"
#include "core/volume.hxx"
#include "core/face.hxx"
#include "core/edge.hxx"
#include "core/parent.hxx"
#include "core/communic.hxx"
#include "core/util.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ctype.h>

#define FNMAX 80 


//__________________________________ Mesh _____________________________________


Mesh :: Mesh () 
    : Vector<Element*>(), GarbagedObject()
{
    // Constructor. Initializes the receiver to an empty mesh.
    // Programming note: C++ is a shit: attribute 'nbCoordinates' below cannot
    //    be set to 'GetDimension()', even if GetDimension() is declared as a
    //    non-pure virtual method.

    skeleton                 = NULL ;

    coordinates              = NULL ;
    jacobian                 = NULL ;
    invJacobian              = NULL ;
    normal                   = NULL ;

    nbCoordinates            = UNINITIALIZED_INTEGER ;
    defaultParentElement     = ParentElementGLL::GetParentElementOfDegree(1) ;
    inhibitedCoordinates     = false ;

    internalCommunicators    = new Vector<Communicator*>() ;
    sendReceiveCommunicators = new Vector<Communicator*>() ;
    receiveCommunicators     = new Vector<Communicator*>() ;
    ComputeLocalElements();
}


Mesh :: ~Mesh () 
{
    // Destructor.

    Element *elem ;
    int     i ;

    for (i=1 ; i<=size ; i++) 
    {
	elem = At(i) ;
	unref(elem) ;
    }

    unref(skeleton) ;

    for (i=1 ; i<=internalCommunicators->GetSize() ; i++)
	delete internalCommunicators->At(i) ;
    delete internalCommunicators ;

    for (i=1 ; i<=sendReceiveCommunicators->GetSize() ; i++)
	delete sendReceiveCommunicators->At(i) ;
    delete sendReceiveCommunicators ;

    for (i=1 ; i<=receiveCommunicators->GetSize() ; i++)
	delete receiveCommunicators->At(i) ;
    delete receiveCommunicators ;

    // does not delete or unref attributes coordinates, jacobian, invJacobian and
    // normal: this is done in method BreakCyclicReferences.
}


void Mesh :: BreakCyclicReferences ()
{
    // If the receiver is not referenced any longer, except in cyclic way by its
    // own basic fields, and if none of its basic fields is referenced any
    // longer, except from the receiver, then deletes these basic fields to allow
    // the destruction of the receiver.
  
    int cyclicReferences ;

    cyclicReferences = 0 ;
    if (coordinates)
	cyclicReferences++ ;
    if (jacobian)
	cyclicReferences++ ;
    if (invJacobian)
	cyclicReferences++ ;
    if (normal)
	cyclicReferences++ ;

    if (counter == cyclicReferences) {
	// only the receiver's basic fields point to the receiver

	// does anyone (beside the receiver) point to one the basic fields ? if
	// yes, do not delete the receiver (nor its basic fields)

	if ((coordinates && coordinates->counter > 1) ||
	    (jacobian    && jacobian   ->counter > 1) || 
	    (invJacobian && invJacobian->counter > 1) ||
	    (normal      && normal     ->counter > 1))
	    return ;

	// initialize temporarily the counter to a large values, to prevent the
	// receiver from being deleted by one of its basic fields
	counter = 1000 ;

	// delete the basic fields
	delete coordinates ;
	delete jacobian ;
	delete invJacobian ;
	delete normal ;

	// set the counter to 0, to delete the receiver 
	counter = 0 ;
    }
}

std::vector<int> const& Mesh :: GetLocalElements() const {
    return localElements;
}

void Mesh :: ComputeLocalElements() {
    localElements.clear();
    for (int i=1 ; i<=size ; i++) {
        if(At(i)->GetProcessorNumber() == Mpi_rank) {
            localElements.push_back(i);
        }
    }
}

Mesh* Mesh :: CreateMesh (int dim)
{
    // Returns a new mesh of dimension 'dim'.

# ifdef REQUIRE
    Require("valid argument 'dim'", dim >= 0 && dim <= 3) ;
    InMethod("Mesh::CreateMesh(dim)") ;
# endif

    if (dim == 0)
	return new Mesh0D() ;

    else if (dim == 1)
	return new Mesh1D() ;

    else if (dim == 2)
	return new Mesh2D() ;

    else
	return new Mesh3D() ;
}

void Mesh :: DeleteBasicFields ()
{
    // Deletes the basic fields (coordinates, jacobian, etc) of the receiver.
    // Useful for modifying their interpolation.

    if (coordinates && coordinates->counter == 1)
	delete coordinates ;
    else
	unref(coordinates) ;
    coordinates = NULL ;

    if (jacobian && jacobian->counter == 1)
	delete jacobian ;
    else
	unref(jacobian) ;
    jacobian = NULL ;

    if (invJacobian && invJacobian->counter == 1)
	delete invJacobian ;
    else
	unref(invJacobian) ;
    invJacobian = NULL ;

    if (normal && normal->counter == 1)
	delete normal ;
    else
	unref(normal) ;
    normal = NULL ;
}


void Mesh :: DeleteInvJacobian ()
{
    // Deletes the inverse-jacobian field of the receiver.
    // Useful for making sure that this field has the same interpolation as the
    // coordinate field.

    if (invJacobian->counter == 1)
	delete invJacobian ;
    else
	unref(invJacobian) ;

    invJacobian = NULL ;
}


void Mesh :: DeleteJacobian ()
{
    // Deletes the jacobian field of the receiver.
    // Useful for making sure that this field has the same interpolation as the
    // coordinate field.

    if (jacobian->counter == 1)
	delete jacobian ;
    else
	unref(jacobian) ;

    jacobian = NULL ;
}


void Mesh :: DispatchElements ()
{
    // Assigns a processor number to every element of the receiver.
    //
    // A non-optimized algorithm is used: roughly said, the first "m" elements
    // are assigned to the first processor, the next "m" elements are assigned to
    // the second processor, etc.
    //
    // This method applies whichever the method used for dispatching the elements
    // of the receiver.
    //
    // On a scalar implementation, dispatching is not required (i.e., this method
    // has no effect). 

    int i, k, m, r, p ;

    m = size / Mpi_nbProcessors ;
    r = size % Mpi_nbProcessors ;

    k = 1 ;
    for (p=0 ; p<Mpi_nbProcessors ; p++) {
	// assign m elements to processor p
	for (i=1 ; i<=m ; i++)
	    At(k++)->SetProcessorNumber(p) ;
	// if necessary, assign 1 more element to processor p
	if (r) {
	    At(k++)->SetProcessorNumber(p) ;
	    r-- ;
	}
    }
    ComputeLocalElements();
}

Mesh* Mesh :: Except (Mesh* m)
{
    // Returns a new mesh containing those of spectral elements of the receiver,
    // except those that also belong to 'm'.

    Mesh    *answer ;
    Element *elem ;
    int     i ;

    answer = CreateMesh(GetDimension()) ;

    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	if (! m->Has(elem))
	    answer->Put(elem) ;
    }

    answer->InterpolateCoordinatesLike(this) ;

    answer -> ComputeLocalElements();
    return answer ;
}


Mesh* Mesh :: Extract (int dir, real c)
{
    // Returns a new mesh containing those of spectral elements of the receiver
    // whose barycenter has its 'dir'-th coordinate equal to 'c'.
    // For example, in a 3D case, it can be used for extracting all elements 
    // lying on the plane Y=c.

# ifdef REQUIRE
    Require ("valid argument 'dir'", dir >= 1 && dir <= 3) ;
    InMethod("Mesh::Extract(dir,coord)") ;
# endif

    Mesh       *answer ;
    Element    *elem ;
    Point      *bary ;
    real       cBary ;
    int        i ;
    const real TOL = 1.e-8 ;

    answer = CreateMesh(GetDimension()) ;

    for (i=1 ; i<=size ; i++) {
	elem  = At(i) ;
	bary  = elem->GetBarycenter() ;
	cBary = bary->GetCoordinate(dir) ;
	if (fabs(cBary-c) < TOL)
	    answer->Put(elem) ;
	unref(bary) ;
    }

    answer->InterpolateCoordinatesLike(this) ;

    answer -> ComputeLocalElements();
    return answer ;
}


Communicator* Mesh :: GetInternalCommunicatorTo (Mesh* m)
{
    // Returns the internal communicator to 'm'. Creates it if the receiver has 
    // not made it yet.

    Communicator *comm ;
    int          i ;

    // search it among the existing communicators
    for (i=1 ; i <= internalCommunicators->GetSize() ; i++) {
	comm = internalCommunicators->At(i) ;
	if (comm->GetMesh2() == m)
	    return comm ;
    }

    // not found, so create it
    comm = new Communicator(this,m,INTERNAL_COMMUNICATOR) ;
    internalCommunicators->Put(comm) ;
    return comm ;
}


Communicator* Mesh :: GetReceiveCommunicatorTo (Mesh* m)
{
    // Returns the receive communicator to 'm'. Creates it if the receiver has 
    // not made it yet.

    Communicator *comm ;
    int          i ;

    // search it among the existing communicators
    for (i=1 ; i <= receiveCommunicators->GetSize() ; i++) {
	comm = receiveCommunicators->At(i) ;
	if (comm->GetMesh2() == m)
	    return comm ;
    }

    // not found, so create it
    comm = new Communicator(this,m,RECEIVE_COMMUNICATOR) ;
    receiveCommunicators->Put(comm) ;
    return comm ;
}


Communicator* Mesh :: GetSendReceiveCommunicatorTo (Mesh* m)
{
    // Returns the send-receive communicator to 'm'. Creates it if the receiver 
    // has not made it yet.

    Communicator *comm ;
    int          i ;

    // search it among the existing communicators
    for (i=1 ; i <= sendReceiveCommunicators->GetSize() ; i++) {
	comm = sendReceiveCommunicators->At(i) ;
	if (comm->GetMesh2() == m)
	    return comm ;
    }

    // not found, so create it
    comm = new Communicator(this,m,SEND_RECEIVE_COMMUNICATOR) ;
    sendReceiveCommunicators->Put(comm) ;
    return comm ;
}


FlatField* Mesh :: GetCoordinates ()
{
    // Returns the coordinate field of the receiver. Constructs it if it does not
    // exist yet.

    if (coordinates == NULL) {
	coordinates = new FlatField(GetNbCoordinates(),this,defaultParentElement) ;
	coordinates->SetToCoordinates() ;
	inhibitedCoordinates = false ;
    }
    else if (inhibitedCoordinates) {
	// the values of the coordinate field have been deleted for saving space
	Printf("Reallocating coordinates of %dD mesh %X\n",GetDimension(),this) ;
	coordinates->ReallocateValues() ;
	coordinates->SetToCoordinates() ;
	inhibitedCoordinates = false ;
    }

    return coordinates ;
}


int Mesh :: GetDimension ()
{
    // Issues an error message.

    NotImplemented("Mesh::GetDimension()") ;
    return 0 ;
}


Element* Mesh :: GetElementCloserTo (Point* p)
{
    // Among the elements of the receiver, returns the one whose barycenter is
    // closer to 'p'.
    // Returns NULL if the receiver is empty.

    Element *elem, *answer ;
    Point   *bary ;
    real    dist, d ;
    int     i ;

    answer = NULL ;
    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	bary = elem->GetBarycenter() ;
	if (answer == NULL) {                      // initialize to first element
	    dist   = bary->GetDistanceTo(p) ;
	    answer = elem ;
	}
	else {                                     // any subsequent element
	    d = bary->GetDistanceTo(p) ; 
	    if (d < dist) {
		dist   = d ;
		answer = elem ;
	    }
	}
	unref(bary) ;
    }
    return answer ;
}


Mesh* Mesh :: GetElementsEncompassing (Element* elem)
{
    // Returns a new mesh containing, among the elements of the receiver, those
    // which include 'elem'.
    // For example, if 'elem' is an edge, returns:
    //  - if the receiver is 1D: the edges identical to 'elem' (same pointer);
    //  - if the receiver is 2D: the faces whose set of edges contains 'elem'.

    Mesh    *answer ;
    Element *elem1 ;
    int     i ;

    answer = CreateMesh(GetDimension()) ;
    for (i=1 ; i<=size ; i++) {
	elem1 = At(i) ;
	if (elem->IsIncludedIn(elem1))
	    answer->Put(elem1) ;
    }
    answer->ComputeLocalElements();
    return answer ;
}
  

FlatField* Mesh :: GetInvJacobian ()
{
    // Returns the matrix of the partial derivatives of the local coordinates
    // (r,s,t) with respect to the cartesian coordinates (X,Y,Z).
    //
    // Returns for a
    //   0D mesh:  an error message
    //
    //   1D mesh:  [ dr/dX ]
    //
    //   2D mesh:  [ dr/dX  dr/dY ]
    //             [ ds/dX  ds/dY ]
    //
    //   3D mesh:  [ dr/dX  dr/dY  dr/dZ ]
    //             [ ds/dX  ds/dY  ds/dZ ]
    //             [ dt/dX  dt/dY  dt/dZ ]
    //
    // The same interpolation as for the coordinate field is used.
    //
    // The components are stored rows by rows. For example, in 2D, the 4
    // components are stored in the following order:
    //   (dr/dX , dr/dY , ds/dX , ds/dY).
    //
    // See the warning in method GetJacobian().

# ifdef REQUIRE
    Require ("mesh of non-nul dimension", GetDimension() > 0) ;
    Require ("nb coordinates = dim", GetNbCoordinates() == GetDimension()) ;
    InMethod("Mesh::GetInverseJacobian()") ;
# endif

    int ncomp ;

    if (invJacobian == NULL) {
	ncomp       = GetDimension() * GetDimension() ;
	invJacobian = GetCoordinates()->DuplicateEmpty(ncomp) ;
	invJacobian->SetToInvJacobian() ;
    }

    return invJacobian ;
}


FlatField* Mesh :: GetJacobian ()
{
    // Returns the jacobian field of the receiver. Creates it if it does not
    // exist yet.
    //
    // This field possesses 1 component. It is initialized to:
    //
    //   . case 1: 1D element in a 1D space
    //       J = | dX/dr |
    //
    //   . case 2: 2D element in a 2D space
    //       J = | dX/dr  dX/ds |
    //           | dY/dr  dY/ds |
    //
    //   . case 3: 3D element in a 3D space
    //       J = | dX/dr  dX/ds  dX/dt |
    //           | dY/dr  dY/ds  dY/dt |
    //           | dZ/dr  dZ/ds  dZ/dt |
    //
    //   . case 4: 1D element in a 2D space
    //       J = norm of (-dY/dr, dX/dr)
    //
    //   . case 5: 2D element in a 3D space
    //       J = norm of the cross product of (dX/dr, dY/dr, dZ/dr)
    //                                    and (dX/ds, dY/ds, dZ/ds)
    //
    // Cases 1, 2 and 3 deal with the standard situations.
    // Case 4 is typical of an element on the skin of a 2D mesh.
    // Case 5 is typical of an element on the skin of a 3D mesh.
    // The same interpolation as for the coordinate field is used.

# ifdef REQUIRE
    Require ("mesh of non-nul dimension", GetDimension() > 0) ;
    Require ("valid type of mesh", GetNbCoordinates() == GetDimension() ||
	     GetNbCoordinates() == GetDimension()+1) ;
    InMethod("Mesh::GetJacobian()") ;
# endif

    if (jacobian == NULL) {
	jacobian = GetCoordinates()->DuplicateEmpty(1) ;
	jacobian->SetToJacobian() ;
    }

    return jacobian ;
}


int Mesh :: GetNbCoordinates ()
{
    // Returns the number of components (1, 2 or 3) of the coordinate field of 
    // the receiver.
    //
    // Observe that the number of (X,Y,Z) coordinates is not necessarily 
    // identical to the dimension of the receiver, i.e., to the number of (r,s,t)
    // paramaters of its elements.

    if (nbCoordinates == UNINITIALIZED_INTEGER)
	nbCoordinates = GetDimension() ;

    return nbCoordinates ;
}


FlatField* Mesh :: GetNormal (RealVector* external)
{
    // Returns the unit normal vector "n" on the contour of the receiver.
    //
    // If the normal does not exist yet, computes it. Makes sure that at every
    // dof "n" has same orientation as 'external', i.e., 
    //   n . external >= 0.
    //
    // This normal is defined for 1D and 2D meshes only.

# ifdef REQUIRE
    Require ("nbCoordinates == dimension+1",
	     nbCoordinates == UNINITIALIZED_INTEGER
	     || nbCoordinates == GetDimension()+1) ;
    InMethod("Mesh::GetNormal()") ;
# endif

    if (normal == NULL) {
	normal = GetCoordinates()->DuplicateEmpty() ;
	normal->SetToNormal() ;
	normal->OrientLike(external) ;
    }

    return normal ;
}


Mesh* Mesh :: GetSkeleton ()
{
    // Returns a new mesh containing the subelements of all elements of the
    // receiver. The dimension of the answer is the dimension of the receiver,
    // minus 1.
    //
    // For every element assigned to a processor p, each of its subelements 
    // (e.g., each of the edges of a face) is assigned to the same processor p,
    // unless it has already be assigned to another processor.

# ifdef REQUIRE
    Require("dimension > 0", GetDimension() > 0) ;
    InMethod("Mesh::GetSkeleton()") ;
# endif

    Element *elem, *subelem ;
    int     i, j, n, p ;

    if (! skeleton) {
	Printf("\nStart computing %dD skeleton of mesh %X\n",
	       GetDimension()-1,this) ; 
	skeleton = CreateMesh(GetDimension()-1) ;
	for (i=1 ; i<=size ; i++) {
	    elem = At(i) ;
	    n    = elem->GetNbSubelements() ;
	    for (j=1 ; j<=n ; j++) {
		subelem = elem->GetSubelement(j) ;
		if (! skeleton->Has(subelem))
		    skeleton->Put(subelem) ;
		p = subelem->GetProcessorNumber() ;
		if (p == UNINITIALIZED_INTEGER)
		    subelem->SetProcessorNumber(elem->GetProcessorNumber()) ;
	    }
	}
	skeleton->InterpolateCoordinatesLike(this) ;
	Printf("OK %dD skeleton\n",skeleton->GetDimension()) ;
    }

    skeleton -> ComputeLocalElements();
    return skeleton ;
}


Mesh* Mesh :: GetSkin ()
{
    // Returns a new mesh containing the subelements on the boundary of the
    // receiver. The dimension of the new mesh is the dimension of the receiver,
    // minus 1.

    Mesh             *answer ;
    Element          *elem, *subelem ;
    Vector<int>      *appearances ;
    Vector<Element*> *subelements ;
    int              i, j, k, n, m, num ;
    const int        MAX_SUB = 20 ;      // max number of subelements per element
  
    m           = size * MAX_SUB ;       // upperbound of total number of subel-
    subelements = new Vector<Element*>(m) ;                           // ements
    appearances = new Vector<int>(m) ;
    subelements->SetValues(NULL) ;
    appearances->SetValues(0) ;

    num = 0 ;
    for (i=1 ; i<=size ; i++) {                    // for every element
	elem = At(i) ;
	n    = elem->GetNbSubelements() ;
	for (j=1 ; j<=n ; j++) {                     // for each of its subelements
	    subelem = elem->GetSubelement(j) ;
	    k       = subelements->GetIndexOf(subelem) ;
	    if (k == 0) {
		k = ++num ;
		subelements->At(k) = subelem ;
	    }
	    appearances->At(k) = appearances->At(k) + 1 ;
	}
    }

    answer = CreateMesh(GetDimension()-1) ;
    for (i=1 ; i<=num ; i++)
	if (appearances->At(i) == 1)
	    answer->Put(subelements->At(i)) ;

    delete subelements ;
    delete appearances ;

    answer->InterpolateCoordinatesLike(this) ;

    answer -> ComputeLocalElements();
    return answer ;
}


void Mesh :: InhibitCoordinates ()
{
    // In order to save memory space, deallocates the values of the coordinate
    // field. 
    //
    // The interpolation of the coordinate field is not deleted, so the values
    // can be retrieved at any time, using the usual method GetCoordinates().
  
    if (coordinates && !inhibitedCoordinates) {
	coordinates->DeallocateValues() ;
	inhibitedCoordinates = true ;
    }
}


void Mesh :: InterpolateCoordinatesLike (Mesh* m)
{
    // Makes sure that the coordinate field (as well as the jacobian field, the 
    // inverse jacobian field and the normal field) of the receiver is
    // interpolated in the same way as the coordinate field of 'm'.
    //
    // The receiver must be included in 'm', and can be of a smaller dimension
    // than 'm'. This method is typically used for making sure that the skin of
    // a mesh is interpolated like 'm', and thus a Dirichlet condition defined on
    // the skin will be correctly imposed on a field defined on 'm'.
    //
    // Note that this also implies that the receiver will have the same number of
    // coordinates as 'm'. For example, if 'm' is a 3D mesh and the receiver is
    // its 2D skin, then the receiver will have 3 coordinates, not the default 2.

    DeleteBasicFields() ;

    nbCoordinates = m->GetNbCoordinates() ;
    coordinates   = new FlatField(nbCoordinates,m->GetCoordinates(),this) ;
    coordinates->SetToCoordinates() ;
}


boolean Mesh :: IsDispatched ()
{
    // Returns 'true' if every element of the receiver has been attached to a
    // processor, else returns 'false'.
    // On a scalar implementation, always returns 'true'.

# ifdef PARALLEL    // parallel implementation

    Element *elem ;
    int     i ;

    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	if (elem->GetProcessorNumber() == UNINITIALIZED_INTEGER)
	    return false ;
    }

    return true ;

# else              // scalar implementation

    return true ;

# endif
}


void Mesh :: Print ()
{
    // Prints the receiver on standard output.

    Element* elem ;
    int      i ;

    printf("Mesh %d of dimension %d with %d elements:\n",
	   this,GetDimension(),size) ;

    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	printf("Element %d:\n",i) ;
	elem->Print() ;
    }
}
 

void Mesh :: PrintBarycenters ()
{
    // A debugging tool.
    // Prints on standard output the barycenter of every element in the receiver.
    // Less verbose than method 'Print'.

    Element *elem ;
    Point   *bary ;
    int     i ;

    if (Mpi_rank == 0) {
	printf("\nBarycenter of every element of %dD mesh %X:\n",
	       GetDimension(),this) ;
	for (i=1 ; i<=size ; i++) {
	    elem = At(i) ;
	    printf ("  elem %2d (on proc. %d): ",i,elem->GetProcessorNumber()) ;
	    bary = elem->GetBarycenter() ;
	    bary->Print() ;
	    unref(bary) ;
	}
    }
}


void Mesh :: PrintDispatching ()
{
    // Prints on standard output the number of elements of the receiver that are
    // attached to every processor.
    // Useful for a parallel implementation.
    //
    // In order to know to which processor every element is attached, invoke
    // rather method PrintBarycenters().

    Element *elem ;
    int     i, n ;

    Printf("\nDispatching of elements of %dD mesh %X:\n",GetDimension(),this) ;

    n = 0 ;
    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	if (elem->GetProcessorNumber() == Mpi_rank)
	    n++ ;
    }

    printf("  number of elements on processor %d: %d\n",Mpi_rank,n) ;
}


void Mesh :: PrintPointers ()
{
    // A debugging tool.
    // Prints on standard output the elements in the receiver. Only pointers are
    // printed. Useful for checking if adjacent elements really share their 
    // common data.

    Element *elem ;
    int     i ;

    printf("Mesh %d of dimension %d: ",this,GetDimension()) ;

    for (i=1 ; i<=size ; i++) {
	elem = At(i) ;
	printf ("  elem %2d:\n",i) ;
	elem->PrintPointers() ;
    }
}


void Mesh :: SetDefaultParentElement (ParentElement* p)
{
    // Sets to 'p' the default parent element for interpolating the 4 basic 
    // fields of the receiver (coordinates, jacobian, inverse jacobian, normal).

    DeleteBasicFields() ;

    defaultParentElement = p ;
}


//__________________________________ Mesh0D ___________________________________


Mesh0D :: Mesh0D () 
    : Mesh()
{
    // Constructor. Initializes the receiver to an empty 0D mesh.
}


Mesh0D :: ~Mesh0D ()
{
    // Destructor.
}


Mesh* Mesh0D :: GetSkeleton ()
{
    // Issues an error message.

    Error("Mesh0D::GetSkeleton","nonsense") ;
    return NULL ;
}


Mesh* Mesh0D :: GetSkin ()
{
    // Issues an error message.

    Error("Mesh0D::GetSkin","nonsense") ;
    return NULL ;
}


void Mesh0D :: Put (Element* elem)
{
    // Adds the vertex 'elem' to the receiver.

# ifdef REQUIRE
    Require("'elem' is a vertex", elem->IsVertex()) ;
    InMethod("Mesh0D::Put(elem)") ;
# endif

    ref(elem) ;
    Vector<Element*>::Put(elem) ;
}


//__________________________________ Mesh1D ___________________________________


Mesh1D :: Mesh1D ()
    : Mesh()
{
    // Constructor. Initializes the receiver to an empty 1D mesh.
} 


Mesh1D :: ~Mesh1D ()
{
    // Destructor.
}


Mesh1D* Mesh1D :: ExtractMeshBetween (Point* pA, Point* pB, Point* pC)
{
    // Returns a new 1D mesh containing those of the edges of the receiver which
    // can link 'pA', 'pB' and 'pC' (in that order).
    //
    // 'pA', 'pB' and 'pC' must coincide geometrically with vertices of the
    // receiver. In other words, 'pA' (and 'pB' and 'pC') cannot be "in the 
    // middle" of an edge.
    //
    // Every edge of the receiver is expected to be connected to exactly two
    // edges.

    Mesh1D  *answer, *edgesA, *edges ;
    Edge    *edge ;
    Vertex  *v ;
    Point   *p1, *p2 ;
    int     i ;
    boolean found_pB ;

    edgesA = GetEdgesWith(pA) ;
    if (!edgesA || edgesA->GetSize() != 2)
	Error("Mesh1D::ExtractMeshBetween(pA,pB,pC)","Unexpected 'edges'") ;

    found_pB = false ;
    for (i=1 ; i<=2 ; i++) {             // two possible branches from pA
	answer = new Mesh1D() ;
	edge   = (Edge*) edgesA->At(i) ;   // pick either of the two edges with pA
	answer->Put(edge) ;
	v  = edge->GetVertexAt(pA) ;
	p2 = edge->GetOtherVertex(v)->GetPoint() ;
	if (p2->IsIdenticalTo(pB))
	    found_pB = true ;

	while (! p2->IsIdenticalTo(pC)) {  // stop when pC is reached
	    p1    = p2 ;
	    edges = GetEdgesWith(p1) ;
	    if (edges->GetSize() != 2)
		Error("Mesh1D::GetSkinBetween","Unexpected size of 'edges'") ;
	    edges->Remove(edge) ;            // do not consider the previous edge
	    edge = (Edge*) edges->At(1) ;
	    answer->Put(edge) ;
	    delete edges ;
	    v  = edge->GetVertexAt(p1) ;
	    p2 = edge->GetOtherVertex(v)->GetPoint() ;
	    if (p2->IsIdenticalTo(pB))
		found_pB = true ;
	    if (p2->IsIdenticalTo(pC)) {
		if (found_pB) {
		    delete edgesA ;                // the correct branch (pA,pC) with pB
		    answer->InterpolateCoordinatesLike(this) ;
                    answer->ComputeLocalElements();
		    return answer ;
		}
		else {
		    delete answer ;                // wrong branch (pA,pC) without pB,
		    answer = NULL ;                // so try the other branch i=2
		    break ;
		}
	    }
	}
    }

    Error("Mesh1D::ExtractMeshBetween(pA,pB,pC)","No branch (pA,pB,pC)") ;
    return NULL ;
}


Mesh1D* Mesh1D :: ExtractMeshBetween (Point* pA, Point* pB)
{

    Mesh1D *answer ;
    answer = new Mesh1D() ;

    for (int i=1 ; i<=size ; i++) {
	Edge* edge = (Edge*)At(i) ;
	if ( edge->HasVertexAt(pA) &&
             edge->HasVertexAt(pB) )
        {
	    answer->Put(edge) ;
            break;
        }
    }
   
    answer -> ComputeLocalElements();
    return answer ;
}

Mesh1D* Mesh1D :: GetEdgesWith (Point* p)
{
    // Returns a new mesh containig those of the receiver's edges which have a
    // vertex coinciding with 'p'.

    Mesh1D *answer ;
    Edge   *edge ;
    int    i ;
   
    answer = new Mesh1D() ;

    for (i=1 ; i<=size ; i++) {
	edge = (Edge*)At(i) ;
	if (edge->HasVertexAt(p)) 
	    answer->Put(edge) ;
    }
   
    answer -> ComputeLocalElements();
    return answer ;
}


void Mesh1D :: Put (Element* elem)
{
    // Adds the edge 'elem' to the receiver.
    // If 'elem' has a vertex 'v2' identical to a vertex 'v1' of an edge of the 
    // receiver, then 'elem' will point to 'v1' rather than to 'v2'.

# ifdef REQUIRE
    Require("'elem' is an edge", elem->IsEdge()) ;
    InMethod("Mesh1D::Put(elem)") ;
# endif

    Edge   *edge1, *edge2 ;
    Vertex *v1, *v2 ;
    int    i, j, k ;

    edge2 = (Edge*)elem ;
    for (i=1 ; i<=size ; i++) {
	edge1 = (Edge*)At(i) ;
	for (j=1 ; j<=2 ; j++) {
	    v1 = edge1->GetVertex(j) ;
	    for (k=1 ; k<=2 ; k++) {
		v2 = edge2->GetVertex(k) ;
		if (v2->IsIdenticalTo(v1))
		    edge2->Substitute(v2,v1) ;
	    }
	}
    }

    ref(edge2) ;
    Vector<Element*>::Put(edge2) ;

    unref(skeleton) ;
    skeleton = NULL ;
}


void Mesh1D :: Put (Mesh1D* mesh2)
{
    // Adds the elements of 'mesh2' to the receiver. Whenever an edge 'e2' of
    // 'mesh2' has a vertex 'v2' identical to a vertex 'e1' of an edge of the
    // receiver, 'e2' will point to 'v1' rather than to 'v2'.
    //
    // The receiver and 'mesh2' can either be adjacent, or overlap, or be dis-
    // joint.
    //
    // The elements in mesh2 will be pointed to by both mesh2 and the receiver.
  
    Edge *edge2 ;
    int  i ;

    for (i=1 ; i<=mesh2->GetSize() ; i++) {
	edge2 = (Edge*) mesh2->At(i) ;
	Put(edge2) ;
    }

    unref(skeleton) ;
    skeleton = NULL ;
}


//_________________________________ Mesh2D ____________________________________


Mesh2D :: Mesh2D ()
    : Mesh()
{
    // Constructor. Initializes the receiver to an empty 2D mesh.
} 


Mesh2D :: ~Mesh2D ()
{
    // Destructor.
}

void Mesh2D :: DispatchElementsByBlocks (int nbElem1, int nbElem2)
{
    // Assigns a processor number to every element of the receiver.
    //
    // Unlike method DispatchElements, this is not a general procedure. The
    // receiver is assumed to be structured, in the sense that its elements are
    // assumed to be generated by a nbElem1*nbElem2 distribution on a face
    // element.
    // The elements are assigned by blocks, in order to minimize inter-processor
    // communication.
    //
    // On a scalar implementation, dispatching is not required (i.e., this method
    // has no effect). 

    Element *elem ;
    real    x ;
    int     *nbElemPerProc1, *nbElemPerProc2 ;
    int     i, j, k, nbProc, nbProc1, nbProc2, p, p1, p2, remnant ;

    // 1. Get the number of blocks (i.e, of processors) in directions 1 and 2

    nbProc = Mpi_nbProcessors ;
    x = sqrt(real(nbProc * nbElem1/nbElem2)) ;
    if (ceil(x)-x > x-floor(x)) {          // nbProc1 = closest integer to x
	nbProc1 = int(ceil(x)) ;
	k = -1 ;
    }
    else {
	nbProc1 = int(floor(x)) ;
	k = 1 ;
    }
    if (nbProc1 == 0)
	nbProc1 = 1 ;

    while (true)
	if (nbProc % nbProc1 == 0) {
	    // acceptable block decomposition
	    nbProc2 = nbProc/nbProc1 ;
	    break ;
	}
	else {
	    // non-integer value of nbProc2: increase or reduce nbProc1. try in the
	    // sequence nbProc1+1, nbProc1-1, nbProc1+2, nbProc1-2, etc (if k=1) or
	    // nbProc1-1, nbProc1+1, nbProc1-2, nbProc1+2, etc (if k=-1)
	    // if kkk
	    nbProc1 += k ;
	    k = -k - (k/abs(k)) ;
	}

    // 2. Get the number of elements to assign to every processor
    //    This number may differ from processor to processor

    nbElemPerProc1 = new int[nbProc1] ;
    remnant        = nbElem1 % nbProc1 ;
    for (p=0 ; p<nbProc1 ; p++) {
	if (p < remnant)
	    // the first 'remnant' processors receive 1 more element than the others
	    nbElemPerProc1[p] = nbElem1/nbProc1 + 1 ;
	else
	    nbElemPerProc1[p] = nbElem1/nbProc1 ;
    }

    nbElemPerProc2 = new int[nbProc2] ;
    remnant        = nbElem2 % nbProc2 ;
    for (p=0 ; p<nbProc2 ; p++)
	if (p < remnant)
	    // the first 'remnant' processors receive 1 more element than the others
	    nbElemPerProc2[p] = nbElem2/nbProc2 + 1 ;
	else
	    nbElemPerProc2[p] = nbElem2/nbProc2 ;

    // 3. Assign the elements to the processors, by blocks

    i  = 0 ;
    j  = 1 ;
    p  = 0 ;
    p1 = 0 ;
    p2 = 0 ;
    for (k=1 ; k<=size ; k++) {
	elem = At(k) ;
	i++ ;
	if (i > nbElemPerProc1[p1]) {          // change block in direction 1
	    p++ ;
	    p1++ ;
	    i = 1 ;
	    if (p1 == nbProc1) {                 // switch to next row
		p -= nbProc1 ;
		p1 = 0 ;
		j++ ;
		if (j > nbElemPerProc2[p2]) {      // change block in direction 2
		    p += nbProc1 ;
		    p2++ ;
		    j = 1 ;
		}
	    }
	}
	elem->SetProcessorNumber(p) ;
    }
    ComputeLocalElements();
    /* Valgrind asked to comment the following two lines
        delete nbElemPerProc1 ;
        delete nbElemPerProc2 ;
    */
}


void Mesh2D :: Put (Element* elem)
{
    // Adds the face 'elem' to the receiver.
    // If 'elem' has an edge 'e2' identical to an edge 'e1' of a face of the 
    // receiver, then 'elem' will point to 'e1' rather than to 'e2' (this may
    // require modifying the orientation of that edge in 'elem').

    Face *face1, *face2 ;
    Edge *edge1, *edge2 ;
    int  i, j, k, n1, n2 ;

# ifdef REQUIRE
    Require("'elem' is a face", elem->IsFace()) ;
    InMethod("Mesh2D::Put(elem)") ;
# endif

    face2 = (Face*)elem ;
    n2    = face2->GetNbEdges() ;
    for (i=1 ; i<=size ; i++) {
	face1 = (Face*) At(i) ;
	n1    = face1->GetNbEdges() ;
	for (j=1 ; j<=n1 ; j++) {
	    edge1 = face1->GetEdge(j) ;
	    for (k=1 ; k<=n2 ; k++) {
		edge2 = face2->GetEdge(k) ;
		if (edge2->IsIdenticalTo(edge1))
		    face2->Substitute(edge2,edge1) ;
	    }
	}
    }

    ref(face2) ;
    Vector<Element*>::Put(face2) ;

    unref(skeleton) ;
    skeleton = NULL ;
}


void Mesh2D :: Put (Mesh2D* mesh2)
{
    // Adds the elements of 'mesh2' to the receiver. Whenever a face 'f2' of 
    // 'mesh2' has an edge 'e2' identical to an edge 'e1' of a face of the
    // receiver, 'f2' will point to 'e1' rather than to 'e2' (this may require
    // modifying the orientation of that edge in 'f2').
    //
    // The receiver and 'mesh2' can either be adjacent, or overlap, or be dis-
    // joint.
    //
    // The elements in 'mesh2' will be attribute of both 'mesh2' and the receiver
    // (i.e., they are not duplicated).
  
    Face *face2 ;
    int  i ;

    for (i=1 ; i<=mesh2->size ; i++) {
	face2 = (Face*) mesh2->At(i) ;
	Put(face2) ;
    }

    unref(skeleton) ;
    skeleton = NULL ;
}



//_________________________________ Mesh3D ____________________________________


Mesh3D :: Mesh3D ()
    : Mesh()
{
    // Constructor. Initializes the receiver to an empty 3D mesh.
} 


void Mesh3D :: DispatchElementsByBlocks (int nbElem1, int nbElem2, int nbElem3)
{
    // Assigns every element of the receiver to a processor. The receiver is
    // assumed to be structured, and its elements are assumed to be generated by
    // an nbElem1*nbElem2*nbElem3 distribution on one generation element.
    //
    // The elements are assigned by blocks, in order to minimize inter-processor
    // communication.
    //
    // As opposed to the 2D case, it is also assumed here that there exists a 
    // partition such that all processors have the same number of elements.
    //
    // On a scalar implementation, dispatching is not required (i.e., this method
    // has no effect). 

# ifdef REQUIRE
    Require("more elements than processors", 
	    nbElem1*nbElem2*nbElem3 >= Mpi_nbProcessors) ;
    InMethod("Mesh3D::DispatchElementsByBlocks(nbElem1,nbElem2,nbElem3)") ;
# endif

    Element      *elem ;
    Vector<int*> *partitions ;
    int          *partition ;
    int          i, j, k, ii, nbProc, nbProc1, nbProc2, nbProc3, nbProc23, 
	nbFaces, minFaces, nb1, nb2, nb3, p, p1, p2 ;

    Printf("\nStart block partitioning and dispatching\n") ;

    // 1. Find all possible combinations

    nbProc     = Mpi_nbProcessors ;
    partitions = new Vector<int*>() ;

    for (nbProc1=1 ; nbProc1 <= max(nbProc,nbElem1) ; nbProc1++)
	if (nbProc % nbProc1 == 0 && nbElem1 % nbProc1 == 0) {
	    nbProc23 = nbProc / nbProc1 ;    // nb processors for directions 2 and 3
	    for (nbProc2=1 ; nbProc2 <= max(nbProc23,nbElem2) ; nbProc2++)
		if (nbProc23 % nbProc2 == 0 && nbElem2 % nbProc2 == 0) {
		    nbProc3 = nbProc23 / nbProc2 ;
		    if (nbElem3 % nbProc3 == 0) {        // admissible partition
			Printf("  admissible partition No. %d: %d x %d x %d procs\n",
			       partitions->GetSize()+1,nbProc1,nbProc2,nbProc3) ;
			partition    = new int[3] ;
			partition[0] = nbProc1 ;
			partition[1] = nbProc2 ;
			partition[2] = nbProc3 ;
			partitions->Put(partition) ;
		    }
		}
	}
    if (partitions->GetSize() == 0)
	Error("Mesh3D::DispatchElementsByBlocks()","found no partition") ;

    // 2. Select the best partition, i.e., that with fewer inter-block faces

    for (i=1 ; i <= partitions->GetSize() ; i++) {
	partition = partitions->At(i) ;
	nbProc1   = partition[0] ;
	nbProc2   = partition[1] ;
	nbProc3   = partition[2] ;
	nb1       = nbElem1/nbProc1 ; // number of elements per proc in direction 1
	nb2       = nbElem2/nbProc2 ;
	nb3       = nbElem3/nbProc3 ;
	nbFaces   = (nbProc1-1)*nbElem2*nbElem3 + (nbProc2-1)*nbElem3*nbElem1 
	    + (nbProc3-1)*nbElem1*nbElem2 ;
	Printf("  number of communication faces of partition %d: %d\n",i,nbFaces) ;

	if (i == 1) {
	    ii = 1 ;
	    minFaces = nbFaces ;
	}
	else 
	    if (nbFaces < minFaces) {
		ii = i ;
		minFaces = nbFaces ;
	    }
    }
    Printf("  -> retained partition: partition No. %d\n",ii) ;
    partition = partitions->At(ii) ;
    nbProc1   = partition[0] ;
    nbProc2   = partition[1] ;
    nbProc3   = partition[2] ;
    nb1       = nbElem1 / nbProc1 ;
    nb2       = nbElem2 / nbProc2 ;
    nb3       = nbElem3 / nbProc3 ;

    // 3. Assign the elements to the processors

    i  = 0 ;
    j  = 1 ;
    k  = 1 ;
    p  = 0 ;
    p1 = 0 ; 
    p2 = 0 ; 
    for (ii=1 ; ii<=size ; ii++) {
	elem = At(ii) ;
	i++ ;
	if (i > nb1) {                      // change block in direction 1 
	    p++ ;
	    p1++ ;
	    i = 1 ;
	    if (p1 == nbProc1) {              // switch to next row
		p -= nbProc1 ;
		p1 = 0 ;
		j++ ;
		if (j > nb2) {                  // change block in direction 2
		    p += nbProc1 ;
		    p2++ ;
		    j = 1 ;
		    if (p2 == nbProc2) {          // switch to next layer
			p -= nbProc1*nbProc2 ;
			p2 = 0 ;
			k++ ;
			if (k > nb3) {              // change block in direction 3
			    p += nbProc1*nbProc2 ;
			    k  = 1 ;
			}
		    }
		}
	    }
	}
	elem->SetProcessorNumber(p) ;
    }

    // 4. Cleaning

    for (i=1 ; i <= partitions->GetSize() ; i++)
	delete partitions->At(i) ;
    delete partitions ;

    Printf("OK block partitioning and dispatching\n") ;
    ComputeLocalElements();
}


void Mesh3D :: Put (Element* elem)
{
    // Adds the volume 'elem' to the receiver.
    // If 'elem' has a face 'f2' identical to a face 'f1' of a volume of the 
    // receiver, then 'elem' will point to 'f1' rather than to 'f2' (this may
    // require modifying the orientation of that face in 'elem').

    Volume *volume1, *volume2 ;
    Face   *face1, *face2 ;
    int    i, j, k, n1, n2 ;

# ifdef REQUIRE
    Require("'elem' is a volume", elem->IsVolume()) ;
    InMethod("Mesh3D::Put(elem)") ;
# endif

    volume2 = (Volume*)elem ;
    n2      = volume2->GetNbFaces() ;
    for (i=1 ; i<=size ; i++) {
	volume1 = (Volume*) At(i) ;
	n1      = volume1->GetNbFaces() ;
	for (j=1 ; j<=n1 ; j++) {
	    face1 = volume1->GetFace(j) ;
	    for (k=1 ; k<=n2 ; k++) {
		face2 = volume2->GetFace(k) ;
		if (face2->IsIdenticalTo(face1))
		    volume2->Substitute(face2,face1) ;
	    }
	}
    }

    ref(volume2) ;
    Vector<Element*>::Put(volume2) ;

    unref(skeleton) ;
    skeleton = NULL ;
}

void Mesh3D :: PutWithoutCheck (Element* elem)
{
    // Adds the volume 'elem' to the receiver.
    // If 'elem' has a face 'f2' identical to a face 'f1' of a volume of the 
    // receiver, then 'elem' will point to 'f1' rather than to 'f2' (this may
    // require modifying the orientation of that face in 'elem').

    Volume *volume1, *volume2 ;
    Face   *face1, *face2 ;
    int    i, j, k, n1, n2 ;

# ifdef REQUIRE
    Require("'elem' is a volume", elem->IsVolume()) ;
    InMethod("Mesh3D::Put(elem)") ;
# endif

    volume2 = (Volume*)elem ;
    ref(volume2) ;
    Vector<Element*>::Put(volume2) ;

    unref(skeleton) ;
    skeleton = NULL ;
}

void Mesh3D :: Put (Mesh3D* mesh2)
{
    // Adds the elements of 'mesh2' to the receiver. Whenever a face 'f2' of 
    // 'mesh2' has an edge 'e2' identical to an edge 'e1' of a face of the
    // receiver, 'f2' will point to 'e1' rather than to 'e2' (this may require
    // modifying the orientation of that edge in 'f2').
    //
    // The receiver and 'mesh2' can either be adjacent, or overlap, or be dis-
    // joint.
    //
    // The elements in 'mesh2' will be attribute of both 'mesh2' and the receiver
    // (i.e., they are not duplicated).
  
    Volume *vol2 ;
    int  i ;

    for (i=1 ; i<=mesh2->size ; i++) {
	vol2 = (Volume *) mesh2->At(i) ;
	Put(vol2) ;
    }

    unref(skeleton) ;
    skeleton = NULL ;
}

