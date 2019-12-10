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

// field.cxx

#include "core/field.hxx"
#include "core/mesh.hxx"
#include "core/element.hxx"
#include "core/face.hxx"
#include "core/edge.hxx"
#include "core/parent.hxx"
#include "core/point.hxx"
#include "core/matrix.hxx"
#include "core/postpro.hxx"
#include "core/communic.hxx"
#include "core/parallel.hxx"
#include "core/solver.hxx"
#include "core/precond.hxx"
#include "core/problem.hxx"
#include "core/timinteg.hxx"
#include "core/bc.hxx"
#include "core/param.hxx"
#include <fstream>
#include <cstring>
static int ione = IONE ;
using namespace std ;

//__________________________________ Field ____________________________________


PolynomialChoice  Field :: polynomialChoice = RICH_PC ;
// Default initialization of a class attribute.


Field :: Field ()
    : GarbagedObject ()
{
    // Constructor.

}


Field :: ~Field ()
{
    // Destructor.

}


void Field :: Assemble ()
{
    // Default implementation: does nothing.
}


void Field :: CopyDofsFromVector (RealVector* x)
{
    // Copies into the receiver, from 'x', the values corresponding to the 
    // degrees of freedom of the receiver.

    int i, n ;

    n = GetNbDofs();
    for (i=1 ; i<=n ; i++)
	GetSetDof(1,i,x->At(i));
}


void Field :: CopyDofsToMatrixColumn (FullMatrix* A, int j)
{
    // Copies into the 'j'-th column of 'A' the values corresponding to the 
    // degrees of freedom.

# ifdef REQUIRE
    Require("valid column number 'j'", j>=1 && j<=A->GetNcol()) ;
    InMethod("Field::CopyDofsToMatrixColumn(A,j)") ;
# endif

    int i, n;

    n= GetNbDofs() ;
    for (i=1 ; i<=n ; i++) 
	GetSetDof(0,i,A->At(i,j)) ;
}


void Field :: CopyDofsToVector (RealVector* x)
{
    // Copies into 'x' the values corresponding to the degrees of freedom.

    int i, n;

    n= GetNbDofs();
    for (i=1 ; i<=n ; i++)
	GetSetDof(0,i,x->At(i));
}


void Field :: Distribute ()
{
    // Default implementation: does nothing.
}


void FlatField :: DivideByWeights ()
{
    // Divide the values of the receiver by the weights of the quadrature rule
    // associated to the discretization.

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->DivideByWeights(ic) ;
    }
}


void Field :: EnforceMortarConstraints ()
{
    // Default implementation: does nothing.
}


Field* Field :: FieldDuplicate ()
{
    // Returns an exact copy of the receiver.
    // Remark. A sounder, polymorphic, design would be to call this method
    //         "Duplicate", to declare it as virtual and to redefine it in the
    //         subclasses FlatField and MortaredField. However, C++ demands that
    //         these redefined versions return a Field*, not a FlatField* or a
    //         MortaredField*; this induces many ugly and potential dangerous
    //         casts.

    FlatField     *flatAnswer ;
    MortaredField *mortaredAnswer ;
  
    if (IsFlatField()) {
	flatAnswer = ((FlatField*)this)->Duplicate() ;
	flatAnswer->CopyFrom(this) ;
	return flatAnswer ;
    }
    else {
	mortaredAnswer = ((MortaredField*)this)->Duplicate() ;
	mortaredAnswer->CopyFrom(this) ;
	return mortaredAnswer ;
    }
}


Field* Field :: FieldDuplicateEmpty (int ncomp)
{
    // Returns a copy of the receiver, with 'ncomp' components, initialized to 0.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.
    //
    // See the remark in method FieldDuplicate().

    FlatField     *flatAnswer ;
    MortaredField *mortaredAnswer ;
  
    if (IsFlatField()) {
	flatAnswer = ((FlatField*)this)->DuplicateEmpty(ncomp) ;
	return flatAnswer ;
    }
    else {
	mortaredAnswer = ((MortaredField*)this)->DuplicateEmpty(ncomp) ;
	return mortaredAnswer ;
    }
}


Field* Field :: FieldGetWork (int ncomp)
{
    // Returns a copy of the receiver, with 'ncomp' components, initialized to 0.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.
    // The values of the answer are initialized to garbage.
    //
    // About the name of this method, see the remark in method FieldDuplicate().

    FlatField     *flatAnswer ;
    MortaredField *mortaredAnswer ;
  
    if (IsFlatField()) {
	flatAnswer = ((FlatField*)this)->GetWork(ncomp) ;
	return flatAnswer ;
    }
    else {
	mortaredAnswer = ((MortaredField*)this)->GetWork(ncomp) ;
	return mortaredAnswer ;
    }
}


PolynomialChoice Field :: GetPolynomialChoice ()
{
    // Returns the value of the class variable 'polynomialChoice'.

    return polynomialChoice ;
}


boolean Field :: IsFlatField ()
{
    // Default: returns 'false', because the receiver is not a flat field.

    return false ;
}


boolean Field :: IsMortaredField ()
{
    // Default: returns 'false', because the receiver is not a mortared field.

    return false ;
}


void Field :: SetName (string s)
{
    // Sets the receiver's name to 's'. Can be useful for debugging and printing.

    name = s ;
}


void Field :: SetPolynomialChoice (PolynomialChoice pc)
{
    // Sets the class attribute 'polynomialChoice' to 'pc', i.e., to RICH_PC or
    // to POOR_PC.
    // Every time an ambiguity is raised about which polynomial expansion must be
    // chosen when creating a mortar, the choice will be governed by 'pc'.

    polynomialChoice = pc ;
}


//_________________________________ FlatField _________________________________


Vector<int>* FlatField :: freeIndices = new Vector<int>() ;
// Initializes a class attribute.


// Initializes some temporary class attributes.
int FlatField :: nbFields0D     = 0 ;
int FlatField :: nbFields1D     = 0 ;
int FlatField :: nbFields2D     = 0 ;
int FlatField :: nbFields3D     = 0 ;
int FlatField :: nbComponents0D = 0 ;
int FlatField :: nbComponents1D = 0 ;
int FlatField :: nbComponents2D = 0 ;
int FlatField :: nbComponents3D = 0 ;


FlatField :: FlatField ()
    : Field()
{
    // Protected constructor, to be used only internally, e.g. in methods for
    // duplicating fields.

    mesh           = NULL ;
    nbComponents   = UNINITIALIZED_INTEGER ;
    index          = UNINITIALIZED_INTEGER ;
    workFields     = new Vector<FlatField*>(0) ;
    availability   = new Vector<boolean>(0) ;
    owner          = NULL ;
    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}


FlatField :: FlatField (int ncomp, Mesh* m, ParentElement* p)
    : Field ()
{
    // Constructor. Initializes the receiver to a field on the elements of 'm', 
    // with 'ncomp' components. Every elementary field has 'p' as parent element
    // in every direction. The values are initialized to 0.

# ifdef REQUIRE
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("Field::Field(ncomp,m,p)") ;
# endif

    Element                *elem ;
    Vector<ParentElement*> *vp ;
    int                    k, n ;

    ref(m) ;
    mesh           = m ;
    nbComponents   = ncomp ;
    index          = GetFreeIndex() ;
    workFields     = new Vector<FlatField*>(0) ;
    availability   = new Vector<boolean>(0) ;
    owner          = NULL ;
    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;

    vp = new Vector<ParentElement*>(mesh->GetDimension()) ;
    vp->SetValues(p) ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->CreateElementaryField(index,ncomp,vp) ;
    }
    delete vp ;

    Count(this) ;
}


FlatField :: FlatField (int ncomp, Mesh* m, Vector<ParentElement*>* vp)
    : Field ()
{
    // Constructor. Initializes the receiver to a field on the elements of 'm',
    // with 'ncomp' components. For every elementary field, the parent element 
    // vp(i) is used as parent element in the i-th direction. The values are 
    // initialized to 0.
    //
    // Remark: every elementary field will point to a copy of 'vp', not to 'vp'
    //         directly.

# ifdef REQUIRE
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("Field::Field(ncomp,m,vp)") ;
# endif

    Element *elem ;
    int     k, n ;

    ref(m) ;
    mesh           = m ;
    nbComponents   = ncomp ;
    index          = GetFreeIndex() ;
    workFields     = new Vector<FlatField*>(0) ;
    availability   = new Vector<boolean>(0) ;
    owner          = NULL ;
    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->CreateElementaryField(index,ncomp,vp) ;
    }

    Count(this) ;
}


FlatField :: FlatField (int ncomp, FlatField* x, Mesh* m)
    : Field ()
{
    // Constructor. Initializes the receiver to a flat field on 'm', with 'ncomp'
    // components. 'm' must be included in the mesh of 'x', and can be of a 
    // smaller dimension than the mesh of 'x'.
    // For every element, the parent elements are chosen as the corresponding 
    // parent elements of the spectral element of the mesh of 'x'.
    // The values of the receiver are initialized to 0.
    //
    // In the case of support meshes of different dimensions, the concept of
    // master and slave spectral elements is used to determine which element is
    // imposing its interpolation rules, when there is more than one candidate
    // element. The choice between the poorest or the richest discretization is
    // governed by the class variable 'mortarChoice' of class MortaredField.
    // For example, if the receiver is defined on a 1D mesh and 'x' on a 2D mesh,
    // an edge of the receiver's mesh can belong to two faces of the mesh of 'x';
    // among these faces, the one with the poorest (resp. richest) interpolation 
    // is chosen for the edge if mortarChoice is set to "poor" (resp. to "rich").
    //
    // This method is used in 2 typical situations:
    // - by the user for creating a new flat field, such as the restriction of
    //   a given field to a smaller-dimensional mesh. Example: the prescribed-
    //   value field of a Dirichlet condition;
    // - by a mortared field when creating its mortar.

# ifdef REQUIRE
    Require("mesh 'm' has correct dimension",
	    m->GetDimension() == x->GetMesh()->GetDimension() ||
	    m->GetDimension() == x->GetMesh()->GetDimension()-1) ;
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("FlatField::FlatField(ncomp,x,m)") ;
# endif

    Element                *elemY, *elemX ;
    Mesh                   *elemsX ;
    Vector<ParentElement*> *vp, *vpMaster ;
    Vector<int>            *degreeMaster ;
    int                    j, k, kk, n, nn, dim, degree ;

    ref(m) ;
    mesh           = m ;
    nbComponents   = ncomp ;
    index          = GetFreeIndex() ;
    workFields     = new Vector<FlatField*>(0) ;
    availability   = new Vector<boolean>(0) ;
    owner          = NULL ;
    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elemY = mesh->At(k) ;

	// select among the elements of mesh(x) those which are identical or
	// encompass elemY. Applies if mesh dimensions are equal or not equal

	elemsX = x->GetMesh()->GetElementsEncompassing(elemY) ;
#   ifdef CHECK
	if (elemsX->GetSize() == 0)
	    Error("FlatField::FlatField(ncomp,x,m)",
		  "mesh of y not included in mesh of x") ;
#   endif
	dim          = mesh->GetDimension() ;
	vp           = NULL ;
	vpMaster     = NULL ;
	degreeMaster = new Vector<int>(dim) ;

	nn = elemsX->GetSize() ;
	for (kk=1 ; kk<=nn ; kk++) {
	    elemX = elemsX->At(kk) ;

	    if (vpMaster == NULL) {                        // initialize (i=1)
		vpMaster = elemY->GetCommonParentElements(elemX,x->GetIndex()) ;
		for (j=1 ; j<=dim ; j++)
		    degreeMaster->At(j) = vpMaster->At(j)->GetPolynomialDegree() ;
	    }
	    else {
		vp = elemY->GetCommonParentElements(elemX,x->GetIndex()) ;
		for (j=1 ; j<=dim ; j++) {                   // loop on every direction
		    degree = vp->At(j)->GetPolynomialDegree() ;
		    if (polynomialChoice==POOR_PC && degree < degreeMaster->At(j)) {
			degreeMaster->At(j) = degree ;
			vpMaster->At(j)     = vp->At(j) ;
		    }
		    else if (polynomialChoice==RICH_PC && degree > degreeMaster->At(j)) {
			degreeMaster->At(j) = degree ;
			vpMaster->At(j)     = vp->At(j) ;
		    }
		}
		delete vp ;
	    }
	}

	elemY->CreateElementaryField(index,ncomp,vpMaster) ;
	delete elemsX ;
	delete degreeMaster ;
	delete vpMaster ;
    }

    Count(this) ;
}
  

FlatField :: ~FlatField ()
{
    // Destructor.

    Element *elem ;
    int     i, k, n ;
    string  s ;
  
    Decount(this) ;

    // 1. elementary fields

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->DeleteElementaryField(index) ;
    }

    // 2. work fields

    for (i=1 ; i<=workFields->GetSize() ; i++) {
	if (availability->At(i) == false) {
	    s = workFields->At(i)->name ;
	    if (!s.empty())
	      Printf("\nWork field %d '%s' is still in use",i,s.c_str()) ;
	    Error("FlatField::~FlatField()","a work field is still in use") ;
	}
	unref(workFields->At(i)) ;            // this should delete the work field
    }
    delete workFields ;
    delete availability ;

    // 3. mesh

    unref(mesh) ;

    freeIndices->At(index) = 0 ; 
}


void FlatField :: Add (Field* x)
{
    // y = y + x .
    // Adds 'x' to the receiver 'y'.
    // HP compiler is a shit: this method should be inherited (and inlined).

    Add(x,ONE) ;
}


void FlatField :: Add (Field* x, real alpha) 
{
    // y = y + alpha x .
    // Adds 'alpha*x' to the receiver 'y'.
  
# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", 
	    nbComponents == x->GetNbComponents()) ;
    InMethod("FlatField::Add(x,alpha)") ;
# endif
  
    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->Add(eX,alpha,ic,ic) ;
    }
}


void FlatField :: Add (Field* x, int icx, int icy)
{
    // y(icy) = y(icy) + x(icx) .
    // Adds the 'icx'-th component of 'x' to the 'icy'-th component of the
    // receiver 'y'.
    // HP compiler is a shit: this method should be inheritable (and inlined).

    Add(x,ONE,icx,icy) ;
}


void FlatField :: Add (Field* x, real alpha, int icx, int icy) 
{
    // y(icy) = y(icy) + alpha x(icx) .
    // Adds the 'icx'-th component of 'x' multiplied by 'alpha' to the 'icy'-th
    // component of the receiver 'y'.

  
# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >=1 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::Add(x,alpha,icx,icy)") ;
# endif
  

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        eY->Add(eX,alpha,icx,icy) ;
    }
}


void FlatField :: Add (real alpha, int icy) 
{
    // y(icy) = y(icy) + alpha .
    // Adds 'alpha' to every value of the 'icy'-th component of the receiver.
    // The case icy=ALL (i.e., icy=0) means: all components are affected.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 0 && icy <= nbComponents) ;
    InMethod("FlatField::Add(alpha,icy)") ;
# endif

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        if (icy == ALL)
            for (ic=1 ; ic<=nbComponents ; ic++)
                eY->Add(alpha,ic) ;
        else
            eY->Add(alpha,icy) ;
    }
}


void FlatField :: AddValuesTo (FlatField* x, NonConform nc)
{
    // Method CopyAddValuesFrom in transposed interpolation mode. It is 
    // implemented in a way that preserves symmetry of operators, i.e., message
    //        y->AddValuesTo(x,nc) ;
    // is equivalent to message
    //        x->AddValuesFrom(y,TRANSPOSED_IM,nc)    (not implemented)
    // meaning that the transpose of the interpolation operator from x to y is 
    // used instead of the interpolation operator from y to x!
    // Copy/AddValuesFrom is thus usually applied to unknown fields whereas
    // Copy/AddValuesTo is intended to be used for operator residual fields.

    x->CopyAddValuesFrom(this,ALL,ALL,ADD_CAM,TRANSPOSED_IM,nc) ;
}

void FlatField :: CopyAddValuesFrom (FlatField* x, int icx, int icy, 
                                     CopyAddMode cam, InterpMode im, 
                                     NonConform nc)
{
    // Manages copy or add values from 'x' to the receiver 'y', with no
    // particular constraints on the meshs of 'x' and 'y', i.e., they can be of 
    // different geometrical dimension or have partial overlap. (This latter
    // feature has not yet been tested).
    //
    // . 'icx' and 'icy' determine which components of 'x' and 'y' are concerned.
    //   If both 'icx' and 'icy' are set to ALL (ie, to 0) the operation involves 
    //   all components of 'x' and 'y'.
    //
    // . 'cam' determines the type of operation:
    //      COPY_CAM: copy x into y,
    //      ADD_CAM:  add x to y;
    //
    // . 'im' determines how to perform interpolation: 
    //      NORMAL_IM:     standard interpolation,
    //      TRANSPOSED_IM: transposed interpolation.
    //   This argument is intended to solve the problem of non equivalence
    //   between the 'inverse' of an interpolation and its transpose, since there
    //   is no equivalence between the two messages:
    //        y->CopyAddValuesFrom(x,...)
    //   and  x->CopyAddValuesTo  (y,...)
    //   The result differs in the choice of matrix interpolation. This non-
    //   equivalence problem does not arise between 0D and 1D fields. See the 
    //   description of method CopyAddValuesTo.
    //
    // . 'nc' determines how to proceed the copy or add operation:
    //      POINTWISE_NC: copy or add using interpolation matrices,
    //      MORTAR_NC:    copy or add using Q operator derived from a weak 
    //                    formulation (mortar). In this case, the interpolation
    //                    matrices mentioned in the hereabove paragraph about
    //                    argument 'im' should be understood as these matrices.
    //   This argument is relevant only in the case of spectral elements of
    //   different dimensions.
    //
    // In the case of elements which overlap partially (eg, a 1D and a 2D
    // element), the copy or add operation is restricted to their intersection.
    //
    // Remarks:
    //
    // - it may happen that data of 'x' should be copied/added into several 
    //   elements of 'y'.
    //
    // - this method applies only to geometrically conformal meshes. A method
    //   like GetElementWithPartialSupport should be implemented for non-
    //   conformal meshes.
    //
    // - in rare cases troubles can be expected: in copy mode, if the mortar 
    //   choice is "rich", if the dimension of the mesh of 'y' is smaller than 
    //   that of the mesh of 'x' (eg, copying faces into an edge), the order may
    //   be relevant. For example, on face1 (resp. face2, edge) the polynomial
    //   degree is 3 (resp., 2, 3); if the copy is performed from face1 (normal)
    //   and then from face2 (reinterpolated) or in the reverse order, the
    //   result is different! So, in methods such as Edge::CopyAddValuesFromFace
    //   copy from a poor element (face2) is not performed. 
    //   However, if both face1 and face2 are of smaller polynomial degree than
    //   edge, then no copy is completed at all! That is why a check is performed
    //   in this method to detect that case. As an example where this case can 
    //   happen, consider a 3D field subject to a bc field on some faces. Since
    //   the copy is deep, edges will copy from faces. The bc field (faces) can
    //   be poorer than the 3D field (edges). The program stops if such case is
    //   encountered.
    //
    // - this method is very heavily used (assembly, distribution, application of
    //   boundary conditions), therefore, for reasons related to numerical
    //   efficiency, the following choice haves been made:
    //   1) the mesh of 'x' and the mesh of 'y' cannot have geometrically 
    //      coincident elements that are not the same element. For example, if 
    //      'y' is 2D and 'x' is 1D, for every face f of the mesh of 'y',
    //      whenever an edge edgeY of f is geometrically identical to an edge
    //      edgeX of 'x', then edgeY and edgeX must be the same object, not two 
    //      objects.
    //      Caveat: this condition is not checked!
    //   2) communicators for improving searches.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 0 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 0 && icy <= nbComponents) ;
    Require("consistent arguments icx and icy",
	    (icx==ALL && icy==ALL) || (icx!=ALL && icy!=ALL)) ;
    Require("same nb components", icx!=ALL || x->nbComponents == nbComponents) ;
    InMethod("Field::CopyAddValuesFrom(x,icx,icy,cam,im,nc)") ;
# endif

    Element      *elemX, *elemY ;
    Element      **pair ;
    Communicator *communicator ;
    int          k, n ; 

    // step 0: initialization
    //         (does nothing in a serial implementation)

    Mpi_ResetCommunication() ;

    // step 1: intra-processor

    communicator = x->mesh->GetInternalCommunicatorTo(mesh) ;
    n  = communicator->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	pair  = communicator->GetPair(k) ;
	elemX = pair[0] ;
	elemY = pair[1] ;
	elemY->CopyAddValuesFrom(elemX,x->index,index,icx,icy,cam,im,nc) ;
    }

    // step 2: inter-processor (send only)
    //         (does nothing in a serial implementation)

    communicator = x->mesh->GetSendReceiveCommunicatorTo(mesh) ;
    n            = communicator->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	pair  = communicator->GetPair(k) ;
	elemX = pair[0] ;
	elemY = pair[1] ;
	elemY->CopyAddValuesFrom(elemX,x->index,index,icx,icy,cam,im,nc,SIZE_RM);
    }

    // step 3: perform communication
    //         (does nothing in a serial implementation)

    Mpi_Communicate() ;

    // step 4: inter-processor (receive only)
    //         (does nothing in a serial implementation)

    communicator = x->mesh->GetReceiveCommunicatorTo(mesh) ;
    n            = communicator->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	pair  = communicator->GetPair(k) ;
	elemX = pair[0] ;
	elemY = pair[1] ;
	elemY->CopyAddValuesFrom(elemX,x->index,index,icx,icy,cam,im,nc,
				 VALUES_RM) ;
    }
}


void FlatField :: CopyFrom (Field* x)
{
    // y = x .
    // Copies 'x' into the receiver 'y'.
  
# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components",
	    nbComponents == x->GetNbComponents()) ;
    InMethod("FlatField::CopyFrom(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->CopyFrom(eX,ic,ic) ;
    }
}


void FlatField :: CopyFrom (Field* x, int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y'.
  
# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    InMethod("FlatField::CopyFrom(x,icx,icy)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        eY->CopyFrom(eX,icx,icy) ;
    }
}


void FlatField :: CopyInterpolateFrom (FlatField* x)
{
    // y = x .
    // Copies 'x' into the receiver 'y'.
    //
    // Although 'y' and 'x' must be defined on the same mesh, their elementary
    // fields need not have the same parent elements: interpolation is provided.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and x have same nb of components", nbComponents==x->nbComponents);
    InMethod("FlatField::CopyInterpolateFrom(x)") ;
# endif

    int ic ;

    for (ic=1 ; ic<=nbComponents ; ic++)
	CopyInterpolateFrom(x,ic,ic) ;
}


void FlatField :: CopyInterpolateFrom (FlatField* x, int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y'.
    //
    // Although 'y' and 'x' must be defined on the same mesh, their elementary
    // fields need not have the same parent elements: interpolation is provided.
  
# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("valid argument 'icx'", icx >=1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::CopyInterpolateFrom(x,icx,icy)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(x->index) ;
        eY->CopyInterpolateFrom(eX,icx,icy) ;
    }
}


void FlatField :: CopyInterpolateTFrom (FlatField* x)
{
    // y = x .
    // Copies 'x' into the receiver 'y'.
    //
    // Although 'y' and 'x' must be defined on the same mesh, their elementary
    // fields need not have the same parent elements: transposed interpolation is
    // provided.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and x have same nb of components", nbComponents==x->nbComponents);
    InMethod("FlatField::CopyInterpolateTFrom(x)") ;
# endif

    int ic ;

    for (ic=1 ; ic<=nbComponents ; ic++)
	CopyInterpolateTFrom(x,ic,ic) ;
}


void FlatField :: CopyInterpolateTFrom (FlatField* x, int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y'.
    //
    // Although 'y' and 'x' must be defined on the same mesh, their elementary
    // fields need not have the same parent elements: transposed interpolation is
    // provided.
  
# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("valid argument 'icx'", icx >=1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::CopyInterpolateTFrom(x,icx,icy)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(x->index) ;
        eY->CopyInterpolateTFrom(eX,icx,icy) ;
    }
}


void FlatField :: CopyValuesFrom (FlatField* x, NonConform nc)
{
    // Method CopyAddValuesFrom in copy mode.

    CopyAddValuesFrom(x,ALL,ALL,COPY_CAM,NORMAL_IM,nc) ;
}


void FlatField :: CopyValuesFrom (FlatField* x, int icx, int icy)
{
    // Method CopyAddValuesFrom in normal interpolation and mortar conformity
    // mode.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >=1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::CopyValuesFrom(x,icx,icy)") ;
# endif 

    CopyAddValuesFrom(x,icx,icy,COPY_CAM,NORMAL_IM,MORTAR_NC) ;
}


void FlatField :: Count (FlatField* x)
{
    // A debugging tool.
    // Takes a record that a new field 'x' has been created.

    int dim ;

    dim = mesh->GetDimension() ;

    if (dim == 0) {
	nbFields0D     += 1 ;
	nbComponents0D += x->nbComponents ;
    }
    else if (dim == 1) {
	nbFields1D     += 1 ;
	nbComponents1D += x->nbComponents ;
    }
    else if (dim == 2) {
	nbFields2D     += 1 ;
	nbComponents2D += x->nbComponents ;
    }
    else {
	nbFields3D     += 1 ;
	nbComponents3D += x->nbComponents ;
    }
}


void FlatField :: Decount (FlatField* x)
{
    // A debugging tool.
    // Takes a record that 'x' is being deleted.

    int dim ;

    dim = mesh->GetDimension() ;

    if (dim == 0) {
	nbFields0D     -= 1 ;
	nbComponents0D -= x->nbComponents ;
    }
    else if (dim == 1) {
	nbFields1D     -= 1 ;
	nbComponents1D -= x->nbComponents ;
    }
    else if (dim == 2) {
	nbFields2D     -= 1 ;
	nbComponents2D -= x->nbComponents ;
    }
    else {
	nbFields3D     -= 1 ;
	nbComponents3D -= x->nbComponents ;
    }
}


void FlatField :: DeallocateValues ()
{
    // Deletes the values of the receiver (but not its interpolation).
    // The purpose is to save memory storage.

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->DeallocateValues() ;
    }
}

real FlatField :: Dot (Field* x)
{
    // dot = y(transposed).x .
    // Returns the scalar product of the receiver 'y' and 'x'.
    // Involves all components of 'y' and 'x'.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("FlatField::Dot(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    real            sum, answer ;
    int             k;

    sum = ZERO ;
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY   = elem->GetElementaryField(index) ;
        eX   = elem->GetElementaryField(((FlatField*)x)->index) ;
        sum += eY->Dot(eX) ;
    }

    answer = Mpi_Sum(sum) ;

    return answer ;
}


FlatField* FlatField :: Duplicate ()
{
    // Returns an exact copy of the receiver.

    FlatField *answer ;

    answer = DuplicateEmpty(ALL) ;
    answer->CopyFrom(this) ;

    return answer ;
}


FlatField* FlatField :: DuplicateEmpty (int ncomp)
{
    // Returns a new field, a copy of the receiver, with 'ncomp' components, and
    // all values initialized to 0.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.

    FlatField *answer ;
    Element   *elem ;
    int       k, n ;

    if (ncomp == ALL)
	ncomp = nbComponents ;

    answer               = new FlatField() ;
    answer->mesh         = mesh ;
    ref(mesh) ;
    answer->nbComponents = ncomp ;
    answer->index        = GetFreeIndex() ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->DuplicateEmptyElementaryField(index,answer->index,ncomp) ;
    }

    Count(answer) ;

    return answer ;
}

void FlatField :: Filter (real alpha, int shift) 
{
    // Applies filter to the receiver according to a method generalized from
    // P. Fischer and J. Mullen's one.
    //
    // F_alpha(u) = alpha I_(Nx-shift).I_(Ny-shift).I_(Nz-shift) u +
    //              (1-alpha) u
    //
    // N+1 is the number of grid points of the mesh corresponding 
    // to the receiver in the corresponding direction. 
    // I_(N-shift)  is the interpolation operator on the grid with only N-shift+1
    // points.
  
# ifdef REQUIRE
    Require(" 0<alpha<1", alpha>=0.0 && alpha<=1.0) ; 
    InMethod("FlatField::Filter(alpha,shift)") ;
# endif

    FlatField* temp ;
    int ic ;
  
    temp = Duplicate();
    temp->ShiftPolynomialDegree(-shift);
    temp->CopyInterpolateFrom(this);
    temp->Multiply(alpha);  
    Multiply(ONE-alpha);
    for (ic=1 ; ic<=GetNbComponents() ; ic++)
	CopyAddValuesFrom(temp,ic,ic,ADD_CAM,NORMAL_IM,POINTWISE_NC);
    
    delete temp;
}
  
real FlatField :: GetEuclidianNorm ()
{
    // Returns the discrete L2 norm of the receiver.
    // As opposed to method GetL2Norm, integration weights are not accounted for
    // here.
  
    return sqrt(Dot(this)) ;
}


int FlatField :: GetFreeIndex ()
{
    // Returns an index that is not yet used by another flat field.
    // Every element uses its index to get the elementary field associated to the
    // receiver.
  
    int answer ;

    answer = freeIndices->GetIndexOf(0) ;  // 0 means: unused
    if (answer == 0) {                     // not found, so enlarge 'freeIndices'
	answer = freeIndices->GetSize() + 1 ;
	freeIndices->Resize(answer) ;
    }
    freeIndices->At(answer) = 1 ;          // 1 means: used

    return answer ;
}


real FlatField :: GetH1Norm ()
{
    // norm = sqrt (integral (y.y) + integral(grad(y) dot grad(y))).
    // Returns the H1 norm of the receiver 'y'.

# ifdef REQUIRE
    Require("not a 0D field", mesh->GetDimension() != 0) ;
    InMethod("FlatField::GetH1Norm()") ;
# endif

    FlatField *gradient ;
    real      normFunction, normGradient ;

    gradient = DuplicateEmpty(nbComponents*mesh->GetDimension()) ;
    gradient->SetToGradient(this) ;
    normFunction = GetL2Norm() ;
    normGradient = gradient->GetL2Norm() ;

    delete gradient ;

    return sqrt(normFunction*normFunction + normGradient*normGradient) ;
}


real FlatField :: GetInteriorInfiniteNorm ()
{
    // Returns the discrete L-infinite norm (max norm) of the receiver. Considers
    // only the interior dofs.
  
    Element *elem ;
    real    answer, norm, elemNorm ;
    int     k;

    norm = ZERO ;
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elemNorm = elem->GetElementaryField(index)->GetInteriorInfiniteNorm() ;
        norm     = max(norm,elemNorm) ;
    }

    answer = Mpi_Max(norm) ;

    return answer ;
}

real FlatField :: GetInfiniteNorm ()
{
    // Returns the discrete L-infinite norm (max norm) of the receiver.
  
    Element *elem ;
    real    answer, norm, elemNorm ;
    int     k;

    norm = ZERO ;
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elemNorm = elem->GetElementaryField(index)->GetInfiniteNorm() ;
        norm     = max(norm,elemNorm) ;
    }

    answer = Mpi_Max(norm) ;

    return answer ;
}


real FlatField :: GetIntegral (int icy)
{
    // Returns the integral of the 'icy'-th component of the receiver, calculated
    // on the receiver's mesh.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::GetIntegral(icy)") ;
# endif

    FlatField *work ;
    real      answer ;


    work = GetWork(1) ;

    work->CopyInterpolateFrom(mesh->GetJacobian()) ;
    work->MultiplyByWeights() ;
    work->Multiply(this,icy,1) ;
    answer = work->SumValues() ;

    Retrieve(work) ;

    return answer ;
}

real FlatField :: GetL2Norm ()
{
    // norm = sqrt (integral (y.y)).
    // Returns the L2 norm of the receiver.
    // Involves all components.
  
# ifdef REQUIRE
    Require("not a 0D field", mesh->GetDimension() != 0) ;
    InMethod("FlatField::GetL2Norm()") ;
# endif

    FlatField *temp, *jacobian ;
    real      norm ;
    int       i ;

    temp = GetWork() ;
    temp->CopyFrom(this) ;

    jacobian = GetWork(1) ;
    jacobian->CopyInterpolateFrom(mesh->GetJacobian()) ;

    for (i=1 ; i<=nbComponents ; i++)
	temp->Multiply(jacobian,1,i) ;
    temp->MultiplyByWeights() ;
    norm = Dot(temp) ;

    Retrieve(temp) ;
    Retrieve(jacobian) ;

    return sqrt(norm) ;
}


int FlatField :: GetMaxNbLocalDofs ()
{
    // Returns the number of dofs (per component) of the elementary field that
    // has more dofs.

    Element *elem ;
    int     k, n, answer, nbDofsElem ;

    answer = 0 ;
    n      = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem       = mesh->At(k) ;
	nbDofsElem = elem->GetElementaryField(index)->GetNbDofs() ;
	answer     = max(answer,nbDofsElem) ;
    }

    return answer ;
}


real FlatField :: GetMeanValue (int icy)
{
    // Returns the mean value of the 'icy'-th component of the receiver.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::GetMeanValue(icy)") ;
# endif

    FlatField *temp ;
    real      integral, area ;

    integral = GetIntegral(icy) ;

    temp = GetWork() ;
    temp->CopyInterpolateFrom(mesh->GetJacobian()) ;
    temp->MultiplyByWeights() ;
    area = temp->SumValues() ;

# ifdef CHECK
    if (area < DEFAULT_TOLERANCE)
	Error("FlatField::GetMeanValue(icy)","area <= 0") ;
# endif

    Retrieve(temp) ;
    return integral/area ;
}


int FlatField :: GetNbDofs ()
{
    // Returns the number of degrees of freedom of the receiver.

    Element *elem ;
    int     k, nbDofsElem, answer ;

    if (nbDofs == UNINITIALIZED_INTEGER) {
	nbDofs = 0 ;
        std::vector<int> const& localElements = mesh->GetLocalElements();
        for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
            k = localElements[iLocal];
            elem = mesh->At(k) ;
            nbDofsElem = elem->GetElementaryField(index)->GetNbDofs()*nbComponents;
            nbDofs    += nbDofsElem ;
	}
    }

    answer = Mpi_Sum(nbDofs) ;

    return answer ;
} 
    

int FlatField :: GetNbInteriorDofs ()
{
    // Returns the number of degrees of freedom of the receiver. The degrees of
    // freedom located on the contour of the elements are not taken into account.

    Element *elem ;
    int     n, k, nbDofsElem ;


    if (nbInteriorDofs == UNINITIALIZED_INTEGER) {
	nbInteriorDofs = 0 ;
	n              = mesh->GetSize() ;
	for (k=1 ; k<=n ; k++) {
	    elem = mesh->At(k) ;
	    nbDofsElem = elem->GetElementaryField(index)->GetNbInteriorDofs()
		* nbComponents ;
	    nbInteriorDofs += nbDofsElem ;
	}
    }

    return nbInteriorDofs ;
} 


FlatField* FlatField :: GetNormByElement (string type)
{
    // Returns a new flat field with one component, with one value per elementary
    // field. This value contains the norm of type 'type' (e.g., H1, L2) of 'x'
    // (all components).
    // This method is expensive (uses probing!).

    FlatField *answer ;

    answer = new FlatField(1,mesh,ParentElementGL::GetParentElementOfDegree(0)) ;
    answer->SetToNormByElement(this,type) ;

    return answer ;
}


FlatField* FlatField :: GetOwner ()
{
    // Returns the owner of the receiver. If the receiver has an owner, every
    // request from the receiver to return a work field will be passed to the
    // owner, in order to limit the proliferation of memory consuming work
    // fields.

    return owner ;
}


void FlatField :: GetSetDof (int iop, int idof, real& val)
{
    // If 'iop' = 0, gets into 'val' the value of the 'idof'-th dof.
    // If 'iop' = 1, sets to 'val' the value of the 'idof'-th dof.

# ifdef REQUIRE
    Require("valid dof number", idof >= 1 && idof <= GetNbDofs()) ;
    InMethod("FlatField::GetSetDof(iop,idof,val)") ;
# endif

    Element *elem ;
    int     k, n, nbDofsElem ;

    n = 0 ;
    for (k=1 ; k <= mesh->GetSize() ; k++) {
	elem       = mesh->At(k) ;
	nbDofsElem = elem->GetElementaryField(index)->GetNbDofs() * nbComponents ;
	if (n+nbDofsElem < idof)      // not the right element yet
	    n += nbDofsElem ;
	else                          // found the right element
	    if (elem->GetProcessorNumber() == Mpi_rank) {
		elem->GetSetDof(index,iop,idof-n,val) ;
		break ;
	    }
    }
}


void FlatField :: GetSetInteriorDof (int iop, int idof, real& val)
{
    // If 'iop' = 0, gets into 'val' the value of the 'idof'-th interior dof.
    // If 'iop' = 1, sets to 'val' the value of the 'idof'-th interior dof.

# ifdef REQUIRE
    Require("valid dof number", idof >= 1 && idof <= GetNbInteriorDofs()) ;
    InMethod("FlatField::GetSetInteriorDof(iop,idof,val)") ;
# endif

    Element *elem ;
    int     k, n, nbDofsElem ;

    n = 0 ;
    for (k=1 ; k <= mesh->GetSize() ; k++) {
	elem       = mesh->At(k) ;
	nbDofsElem = elem->GetElementaryField(index)->GetNbInteriorDofs() 
	    * nbComponents ;
	if (n+nbDofsElem < idof)     // not the right element yet
	    n += nbDofsElem ;
	else                         // found the right element
	    if (elem->GetProcessorNumber() == Mpi_rank) {
		elem->GetSetInteriorDof(index,iop,idof-n,val) ;
		break ;
	    }
    }
}

FlatField* FlatField :: GetWork (int ncomp)
{
    // Returns a copy of the receiver, with 'ncomp' components.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.
    // The values of the answer contain garbage.
    //
    // If the receiver is itself a work field (of its "owner"), then the returned
    // field will be a work field of the owner, not of the receiver. The purpose
    // is to avoid the memory-consuming proliferation of work fields.

# ifdef REQUIRE
    Require("valid argument 'ncomp'", ncomp > 0 || ncomp == ALL) ;
    InMethod("FlatField::GetWork(ncomp)") ;
# endif

    FlatField *w ;
    int       i, n ;

    if (owner)
	w = owner->GetWork(ncomp) ;

    else {
	if (ncomp == ALL)
	    ncomp = nbComponents ;

	n = workFields->GetSize() ;
	for (i=1 ; i<=n ; i++)
	    if (availability->At(i) == true) {
		w = workFields->At(i) ;
		if (w->nbComponents == ncomp) {
		    availability->At(i) = false ;
		    ref(w) ;                                 // now counter should be 2
		    return w ;
		}
	    }
  
	// all work fields are used, so create one
	w = DuplicateEmpty(ncomp) ;
	w->SetOwner(this) ;
	workFields->Put(w) ;
	availability->Put(false) ;
	ref(w) ;                                       // counter = 2
    }

# ifdef CHECK
    if (! w->HasSameInterpolationAs(this))
	Error("FlatField::GetWork(ncomp)","wrong interpolation") ;
# endif

    return w ;
}


boolean FlatField :: HasSameInterpolationAs (Field* x)
{
    // Returns 'true' if:
    //  - the receiver is defined on the same mesh as 'x', and
    //  - on every element the receiver has the same parent elements as 'x'.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    InMethod("FlatField::HasSameInterpolationAs(x)") ;
# endif

    FlatField       *X ;
    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, n ;

    X = (FlatField*)x ;

    if (mesh != X->mesh)
	return false ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	eY   = elem->GetElementaryField(index) ;
	eX   = elem->GetElementaryField(X->index) ;
	if (! eY->HasSameInterpolationAs(eX))
	    return false ;
    }

    return true ;
}


boolean FlatField :: HasSameTypeAs (Field* x)
{
    // Returns 'true' if 'x' is a flat field, else returns 'false'.

    return x->IsFlatField() ;
}


real FlatField :: InteriorDot (Field* x)
{
    // dot = y(transposed).x .
    // Returns the scalar product of the receiver 'y' and 'x'.
    // Involves all components of 'y' and 'x', but only the values at the 
    // interior; for example, for a 2D field, the values on the edges are not
    // taken into account.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("FlatField::InteriorDot(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    real            sum, answer ;
    int             k;

    sum = ZERO ;
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY   = elem->GetElementaryField(index) ;
        eX   = elem->GetElementaryField(((FlatField*)x)->index) ;
        sum += eY->InteriorDot(eX) ;
    }

    answer = Mpi_Sum(sum) ;

    return answer ;
}


void FlatField :: Inverse ()
{
    // y = 1 / y .
    // Substitutes every value of the receiver by its inverse.

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->Inverse() ;
    }
}


boolean FlatField :: IsEqual (FlatField* x)
{
    // Returns 'true' if the receiver and 'x' have identical values, else returns
    // 'false'.

    FlatField *diff ;
    real      norm, normDiff ;

    if (mesh != x->mesh || nbComponents != x->nbComponents 
	|| ! HasSameInterpolationAs(x))
	return false ;

    norm = GetInfiniteNorm() ;
    if (norm < DEFAULT_TOLERANCE) 
	return (x->GetInfiniteNorm() < TWO*DEFAULT_TOLERANCE) ;
    else {
	diff = Duplicate() ;
	diff->Subtract(x) ;
	normDiff = diff->GetInfiniteNorm() ;
	delete diff ;
	return (normDiff/norm < DEFAULT_TOLERANCE) ;
    }
}


boolean FlatField :: IsFlatField ()
{
    // Returns 'true', because the receiver is a flat field.

    return true ;
}

void FlatField :: Multiply (real alpha)
{
    // y = alpha * y.
    // Multiplies every value of the receiver 'y' by 'alpha'.

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic ;

    if (alpha == ONE)
	return ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->Multiply(alpha,ic) ;
    }
}

void FlatField :: Multiply (Field* x)
{
    // y = y * x .
    // Multiplies every value of the receiver 'y' by the vis-a-vis values in 'x'.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have the same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", 
	    nbComponents == x->GetNbComponents()) ;
    InMethod("FlatField::Multiply(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->Multiply(eX,ic,ic) ;
    }
}


void FlatField :: Multiply (Field* x, int icx, int icy)
{
    // y = y * x .
    // Multiplies every value of the 'icy'-th component of the receiver 'y' by 
    // its vis-a-vis in the 'icx'-th component of 'x'.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have the same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >=1 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::Multiply(x,icx,icy)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        eY->Multiply(eX,icx,icy) ;
    }
}


void FlatField :: MultiplyAndAdd (real alpha, Field* x, real beta)
{
    // y = alpha*y + beta*x .
    // Multiplies the receiver 'y' by 'alpha' then adds 'beta * x'.

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", 
	    nbComponents == x->GetNbComponents()) ;
    InMethod("FlatField::MultiplyAndAdd(alpha,x,beta)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(((FlatField*)x)->index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->MultiplyAndAdd(alpha,eX,ic,ic,beta) ;
    }
}


void FlatField :: MultiplyByWeights ()
{
    // Multiplies the values of the receiver by the weights of the quadrature
    // rule associated to the discretization.

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->MultiplyByWeights(ic) ;
    }
}


void FlatField :: OrientLike (RealVector* x)
{
    // Switches the sign of the values of the receiver 'y' wherever the 
    // condition
    //   y.x >= 0 
    // is not met.

# ifdef REQUIRE
    Require("valid dimensions", nbComponents == x->GetSize()) ;
    InMethod("FlatField::OrientLike(x)") ;
# endif

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->OrientLike(x) ;
    }
}


void FlatField :: Print ()
{
    // Prints the receiver on standard output. All values are printed.
    // HP compiler is a shit: this method should be inherited.

    Print(true) ;
}


void FlatField :: Print (boolean showValues)
{
    // Prints the receiver on standard output. Prints the values only if
    // 'showValues' is set to 'true'.

    Element *elem ;
    int     k, n ;

    if (Mpi_rank == 0) {
      if (!name.empty())
	printf("Flat field '%s' %d ",name.c_str(),index) ;
	else
	    printf("Flat field %d ",index) ;

	if (nbComponents <= 1)
	    printf("with %d component, ",nbComponents) ;
	else
	    printf("with %d components, ",nbComponents) ;
	printf("defined on %dD mesh %X\n",mesh->GetDimension(),mesh) ;
    }

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	if (Mpi_rank == elem->GetProcessorNumber())
	    printf("Element %d: ",k) ;
	elem->GetElementaryField(index)->Print(showValues) ;
    }

    if (Mpi_rank == 0)
	printf("\n") ;
}


void FlatField :: Print (Element* elem)
{
    // Prints on standard output the values associated with 'elem'.

# ifdef REQUIRE
    Require("'elem' belongs to the mesh", mesh->Has(elem)) ;
    InMethod("FlatField::Print(elem)") ;
# endif

    elem->GetElementaryField(index)->Print(true) ;
}


void FlatField :: PrintComponent (int icy)
{
    // Prints on standard output the 'icy'-th component of the receiver.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 0 && icy <= nbComponents) ;
    InMethod("FlatField::Print(icy)") ;
# endif
  
    Element *elem ;
    int     k ;

    if (Mpi_rank == 0)
      if (!name.empty())
	    printf("Flat field '%s' with %d components, defined on %dD mesh %X\n",
		   name.c_str(),nbComponents,mesh->GetDimension(),mesh) ;
	else
	    printf("Flat field with %d components, defined on %dD mesh %X\n",
		   nbComponents,mesh->GetDimension(),mesh) ;

    for (k=1 ; k <= mesh->GetSize() ; k++) {
	elem = mesh->At(k) ;
	if (Mpi_rank == 0)
	    printf("Element %d: ",k) ;
	elem->GetElementaryField(index)->PrintComponent(icy) ;
    }
}


void FlatField :: PrintNbFields ()
{
    // A debugging tool.
    // Prints on standard output the number of flat fields currently allocated.

    Printf("\n") ;

    Printf("number of 0D fields:%3d (nb components: %d)\n",
	   nbFields0D,nbComponents0D) ;
    Printf("number of 1D fields:%3d (nb components: %d)\n",
	   nbFields1D,nbComponents1D) ;
    Printf("number of 2D fields:%3d (nb components: %d)\n",
	   nbFields2D,nbComponents2D) ;
    Printf("number of 3D fields:%3d (nb components: %d)\n",
	   nbFields3D,nbComponents3D) ;
}


void FlatField :: PrintNbInstances ()
{
    // A debugging tool. A class method. Prints on standard output the number of
    // created flat fields. 
    // Note: a mortared field on a n-dimensional mesh is made of n+1 flat fields.

    // Obsolete. PrintNbFields() is better.

    int i, nb ;

    nb = 0 ;
    for (i=1 ; i<=freeIndices->GetSize() ; i++)
	if (freeIndices->At(i) == 1)
	    nb++ ;

    if (Mpi_rank == 0)
	printf("\nNumber of flat fields in use: %d\n\n",nb) ;
}


void FlatField :: PrintWork () 
{
    // A debugging tool. Prints on standard output information on the work fields
    // of the receiver.
    // Not parallel yet.

    int i, n ;

    printf("\nWork fields of flat field %d, defined on %dD mesh %X:\n",
           index,mesh->GetDimension(),mesh) ;

    n = workFields->GetSize() ;

    if (n == 0)
	printf("   no work fields\n") ;
    else
	for (i=1 ; i<=n ; i++) {
	    if (availability->At(i) == true)
		printf("%4d  free  ncomp: %d\n",i,workFields->At(i)->nbComponents) ;
	    else
		printf("%4d  used  ncomp: %d\n",i,workFields->At(i)->nbComponents) ;
	}

    printf("\n") ;
}

void FlatField :: ReadFromProc0 (string fName, int icy)
{
    Element    *elem ;
    RealVector *donnees ;
    RealVector *comp ;
    FILE       *file ;
    string      fileName ;
    int         taille=0, k, n ;
    
    fileName = fName ;
    if (Mpi_rank == 0){
      printf(" fileName = %s \n",fileName.c_str()) ;
      file = fopen(fileName.c_str(),"r") ;
	if (!file) {
	    if (Mpi_rank==0) printf("\n Probleme dans ReadFromProc0: avez-vous fait racc00 ?\n");
	    Exit();
	}
    }
    donnees = new RealVector(1) ;
  
    for (k=1 ; k <= mesh->GetSize() ; k++) {
	elem = mesh->At(k) ; 
    
	if (elem->GetProcessorNumber() == Mpi_rank && Mpi_rank != 0){
	    taille = elem->GetElementaryField(index)->GetComponent(icy)->GetSize() ;
	    Mpi_SendInt(&taille,0);
	}
	if (Mpi_rank == 0 && elem->GetProcessorNumber() != 0){
	    Mpi_ReceiveInt(&taille,elem->GetProcessorNumber());
	}
    
	if (Mpi_rank == 0 && elem->GetProcessorNumber() == 0){
	    taille = elem->GetElementaryField(index)->GetComponent(icy)->GetSize() ;
	}
    
	// fprintf(stdout,"elem(proc %d): %d, proc: %d taille: %d\n",elem->GetProcessorNumber(),k,Mpi_rank,taille);
    
	if (Mpi_rank == 0){
	    if (taille != donnees->GetSize())
		donnees->Resize(taille);
      
	    n = fread(donnees->GetValues(),sizeof(real),taille,file) ;
      
	    if (elem->GetProcessorNumber()!=0) 
		Mpi_Send(taille,donnees->GetValues(),1,elem->GetProcessorNumber());
	    else{ 
		comp = elem->GetElementaryField(index)->GetComponent(icy) ;
		comp->CopyFrom(donnees);
	    }
	} 
	else if (elem->GetProcessorNumber() == Mpi_rank && Mpi_rank != 0){
	    comp = elem->GetElementaryField(index)->GetComponent(icy) ;
	    Mpi_Receive(taille,comp->GetValues(),1,0);
	}
	Mpi_Barrier() ;
    }
  
    if (Mpi_rank == 0) fclose(file) ;
    delete donnees ;
}

void FlatField :: PrintInProc0 (string fName, int icy)
{
    Element    *elem ;
    RealVector *donnees ;
    RealVector *comp ;
    FILE       *file ;
    string      fileName ;
    int         size=0, k, n ;
  
    fileName = fName ;
  
    if (Mpi_rank == 0){
      file = fopen(fileName.c_str(),"w") ;
	if (!file) {
	    if (Mpi_rank==0) printf("\n Probleme dans PrintInProc0.\n");
	    Exit();
	}
    }
  
    donnees = new RealVector(1) ;
    donnees->SetValues(0.);
  
    for (k=1 ; k<=mesh->GetSize() ; k++) {
	elem = mesh->At(k) ; 
	
	if (elem->GetProcessorNumber() == Mpi_rank && Mpi_rank != 0){
	    size = elem->GetElementaryField(index)->GetComponent(icy)->GetSize() ;
	    Mpi_SendInt(&size,0);
	}
	if (Mpi_rank == 0 && elem->GetProcessorNumber() != 0){
	    Mpi_ReceiveInt(&size,elem->GetProcessorNumber());
	}
	if (Mpi_rank == 0 && elem->GetProcessorNumber() == 0){
	    size = elem->GetElementaryField(index)->GetComponent(icy)->GetSize() ;
	}
    
	//fprintf(stdout,"elem(proc %d): %d, proc: %d size: %d\n",elem->GetProcessorNumber(),k,Mpi_rank,size);
    
	if (Mpi_rank == 0){
	    if (size != donnees->GetSize())
		donnees->Resize(size);
      
	    if (elem->GetProcessorNumber()!=0) 
		Mpi_Receive(size,donnees->GetValues(),1,elem->GetProcessorNumber());
	    else{ 
		comp = elem->GetElementaryField(index)->GetComponent(icy) ;
		donnees->CopyFrom(comp);
	    }
	    n = fwrite(donnees->GetValues(),sizeof(real),size,file) ;
	} 
	else if (elem->GetProcessorNumber() == Mpi_rank && Mpi_rank != 0){
	    comp = elem->GetElementaryField(index)->GetComponent(icy) ;
	    Mpi_Send(size,comp->GetValues(),1,0);
	}
	Mpi_Barrier() ;
    }
  
    if (Mpi_rank == 0) fclose(file) ;
    delete donnees ;
}

void FlatField :: ReallocateValues ()
{
    // Reallocates the values of the receiver.
    // This method is intended to be invoked for a field whose values have been
    // deallocated (see method DeallocateValues) in order to save memory storage.
    // The values are initialized to 0.

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->AllocateValues() ;
    }
}


void FlatField :: Retrieve (Field* w)
{
    // Makes 'w' available again as a work field identical to the receiver.
    // ('w' may have a different number of components from the receiver).
    // With respect to reference counting, invoking this method is tantamount
    // to calling "unref(w)".

# ifdef REQUIRE
    Require("'w' is a flat field", w->IsFlatField()) ;
    Require("'w' was obtained from 'y'", owner != NULL ||
	    workFields->Has((FlatField*)w)) ;
    Require("same interpolation", HasSameInterpolationAs(w)) ;
    InMethod("FlatField::Retrieve(w)") ;
# endif

    int i, n ;

    if (owner)
	owner->Retrieve(w) ;

    else {
	unref(w) ;                                      // counter should become 1
	n = workFields->GetSize() ;
	for (i=1 ; i<=n ; i++)                          // look for w in the list
	    if (workFields->At(i) == w) {                 // same pointer
#       ifdef CHECK
		if (availability->At(i) == true)
		    Error("FlatField::Retrieve(w)","'w' should be in use") ;
#       endif
		availability->At(i) = true ;
		return ;
	    }
    }
}


void FlatField :: SetLocalDofs (int idof, real val)
{
    // Sets the value of the 'idof'-th local dof of every elementary field of the
    // receiver to 'val'.
    // All components are affected.
    // Does nothing for every elementary field which does not have such 'idof'-th
    // dof.

# ifdef REQUIRE
    Require("valid dof number", idof >= 1 && idof <= GetNbDofs()) ;
    InMethod("FlatField::SetLocalDofs(idof,val)") ;
# endif

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic, nbDofsElem, ncomp ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY         = elem->GetElementaryField(index) ;
        nbDofsElem = eY->GetNbDofs() ;
        if (idof <= nbDofsElem) {
            ncomp = eY->GetNbComponents() ;
            for (ic=1 ; ic<=ncomp ; ic++) 
                eY->GetComponent(ic)->At(idof) = val ;
	}
    }
}


void FlatField :: SetLocalDofs (int idof, FlatField* x)
{
    // Sets the value of the 'idof'-th local dof of every elementary field of the
    // receiver to the corresponding value in 'x'.
    // All components are affected.
    // Does nothing for every elementary field which does not have such 'idof'-th
    // dof.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", x->nbComponents == nbComponents) ;
    InMethod("FlatField::SetLocalDofs(idof,x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, ic, nbDofsElem, ncomp ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY         = elem->GetElementaryField(index) ;
        eX         = elem->GetElementaryField(x->index) ;
        nbDofsElem = eY->GetNbDofs() ;
        if (idof <= nbDofsElem) {
            ncomp = eY->GetNbComponents() ;
            for (ic=1 ; ic<=ncomp ; ic++) 
                eY->GetComponent(ic)->At(idof) = eX->GetComponent(ic)->At(idof) ;
        }
    }
}


void FlatField :: SetOwner (FlatField* x)
{
    // Lets the receiver know that it is a work field of 'x'. Therefore the
    // receiver will never have any work field; instead, if requested a work 
    // field, it will pass the request to 'x'.

# ifdef REQUIRE
    Require("has no owner yet", owner == NULL) ;
    Require("same interpolation", HasSameInterpolationAs(x)) ;
    InMethod("FlatField::SetOwner(x)") ;
# endif

    owner = x ;
}


void FlatField :: SetParentElements (ParentElement* p)
{
    // Changes the discretization rule in every elementary field and in every
    // direction to 'p'.
    // The current values are lost and reinitialized to 0.

    Element *elem ;
    int     k, n, dir ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	for (dir=1 ; dir<=mesh->GetDimension() ; dir++)
	    elem->GetElementaryField(index)->SetParentElement(p,dir) ;
    }
}


void FlatField :: SetParentElements (ParentElement* p, int dir)
{
    // Changes the discretization rule in every elementary field to 'p', in the
    // 'dir'-th direction.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("valid argument 'dir'", dir >= 1 && dir <= mesh->GetDimension()) ;
    InMethod("FlatField::SetParentElements(p,dir)") ;
# endif

    Element *elem ;
    int     k, n ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->GetElementaryField(index)->SetParentElement(p,dir) ;
    }

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}

void FlatField :: SetParentElement (ParentElement* p, int dir, Element* elem)
{
    // Changes the discretization rule in the elementary field of 'elem', in the
    // 'dir'-th direction, to 'p'.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("'elem' belongs to the mesh", mesh->Has(elem)) ;
    Require("valid argument 'dir'", dir >= 1 && dir <= mesh->GetDimension()) ;
    InMethod("FlatField::SetParentElement(p,dir,elem)") ;
# endif

    elem->GetElementaryField(index)->SetParentElement(p,dir) ;

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}


void FlatField :: SetPolynomialDegree (int n)
{
    // Sets to 'n' the polynomial degree of all elementary fields of the 
    // receiver, in all directions.
    //
    // The type of the collocation rule (GL, GLL, etc) need not be the same in
    // all elementary fields, nor in all directions within every elementary 
    // field.
    //
    // The current values are lost and reinitialized to 0.

    Element *elem ;
    int     k, m, dir ;

    m = mesh->GetSize() ;
    for (k=1 ; k<=m ; k++) {
	elem = mesh->At(k) ;
	for (dir=1 ; dir<=mesh->GetDimension() ; dir++)
	    elem->GetElementaryField(index)->SetPolynomialDegree(n,dir) ;
    }

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}


void FlatField :: SetPolynomialDegree (int n, int dir)
{
    // Sets to 'n' the polynomial degree of all elementary fields of the 
    // receiver, in the 'dir'-th directions.
    //
    // The type of the collocation rule (GL, GLL, etc) need not be the same in
    // all elementary fields, nor in all directions within every elementary 
    // field.
    //
    // The current values are lost and reinitialized to 0.

    Element *elem ;
    int     k, m ;

    m = mesh->GetSize() ;
    for (k=1 ; k<=m ; k++) {
	elem = mesh->At(k) ;
	elem->GetElementaryField(index)->SetPolynomialDegree(n,dir) ;
    }

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}


void FlatField :: SetPolynomialDegree (int n, int dir, Element* elem)
{
    // Sets to 'n' the polynomial degree of the elementary field of 'elem', in
    // all directions.
    //
    // The current values of the elementary field of 'elem' are lost and
    // reinitialized to 0.

# ifdef REQUIRE
    Require("'elem' belongs to the mesh", mesh->Has(elem)) ;
    InMethod("FlatField::SetPolynomialDegree(n,dir,elem)") ;
# endif

    elem->GetElementaryField(index)->SetPolynomialDegree(n,dir) ;

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}


void FlatField :: SetTo (Function* f, real t)
{
    // Initializes the values of the receiver to the values of 'f'. A unknown
    // at coordinates (x,y,z) is initialized to f(x,y,z,t).

    int i ;

    for (i=1 ; i<=nbComponents ; i++)
	SetTo(f,i,i,t) ;
}


void FlatField :: SetTo (Function* f, int ic, int icf, real t)
{
    // Initializes the values of the 'ic'-th component receiver to the values of
    // the 'icf'-th component of 'f'.
    // An unknown at coordinates (x,y,z) is initialized to f(x,y,z,t).
    // If the receiver and the mesh's coordinate field have different 
    // interpolation rules, a coordinate field with the same interpolation as the
    // receiver is used, for the sake of accuracy.

    FlatField *coord ;
    Element   *elem ;
    int       k;
    boolean   fast ;

    // 1. Find a coordinate field

    if (HasSameInterpolationAs(mesh->GetCoordinates())) {
	// fast case: use the coordinate field
	fast  = true ;
	coord = mesh->GetCoordinates() ;
    }
    else {
	// slow case: recompute a coordinate field
	fast  = false ;
	coord = GetWork(mesh->GetNbCoordinates()) ;
	coord->SetToCoordinates() ;
    }

    // 2. Initialize the values

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->SetTo(index,coord->index,f,ic,icf,t) ;
    }

    if (! fast)
	Retrieve(coord) ;
}


void FlatField :: SetToAbsoluteValue ()
{
    // Initializes every value of the receiver to its absolute value.

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->SetToAbsoluteValue() ;
    }
}


void FlatField :: SetToCoordinates ()
{
    // Initializes the values of the receiver to the coordinates of the every
    // element.
    // As many components as possible are initialized: for example, if the
    // receiver has 2 components, they will be initialiazed to X and Y,
    // respectively, regardless of the mesh dimension.
    // This method is intended to be used by the coordinate field of the meshes,
    // for initializing its values.

# ifdef REQUIRE
    Require("Valid number of components", nbComponents <= 3) ;
    InMethod("FlatField::SetToCoordinates") ;
# endif

    Element *elem ;
    int     k, n ;

    // Attention: mesh->GetLocalElements() cannot be used here to shorten the
    //   loop. It is unclear why (J. Latt, Jan 2009)
    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	if (elem->GetProcessorNumber() == Mpi_rank) {
            elem->SetToCoordinates(index) ;
        }
    }
}

void FlatField :: SetToConvection (FlatField* x1, FlatField* x2)
{
    // y = (x1 grad x2).
    // Initializes the receiver 'y' to the convection of 'x2' by 'x1',

# ifdef REQUIRE
    Require("y and x1 have same mesh", mesh == x1->mesh) ;
    Require("y and x2 have same mesh", mesh == x2->mesh) ;
    Require("y and x1 are not the same object", this != x1) ;
    Require("y and x2 are not the same object", this != x2) ;
    Require("x2 and y have same nb components",x2->nbComponents == nbComponents);
    InMethod("FlatField::SetToConvection(x1,x2)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eY1, *eY2, *eX1, *eX2, *eInvJacobian ;
    int             i, k, icomp, indexInvJacobian ;

    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX1          = elem->GetElementaryField(x1->index) ;
        eX2          = elem->GetElementaryField(x2->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eY2          = eY->GetWork() ;
        eY1          = eY->GetWork() ;

        eY->SetValues(ZERO) ;
        for (icomp=1 ; icomp<=nbComponents ; icomp++)
            for (i=1 ; i<=mesh->GetDimension() ; i++) {
                eY2->SetToGradient(eX2,icomp,1,i,eInvJacobian) ;
                eY1->CopyInterpolateFrom(eX1,i,1) ;
                eY2->Multiply(eY1) ;
                eY->Add(eY2,1,icomp) ;
            }
        eY->Retrieve(eY2) ;
        eY->Retrieve(eY1) ;
    }
}

void FlatField :: SetToDivergence (FlatField* x, boolean vector)
{
    // y_i = dx_ij/dXj (x tensor)
    // y = dx_j/dXj (x vector)
    // Initializes the receiver 'y' to the divergence (with respect to the 
    // coordinates (X,Y,Z)) of 'x'.
    // All components are concerned.
    // If the parent elements of 'y' differ from those of 'x', re-interpolation
    // from 'x' is provided.
    // The receiver is expected to have enough components for accomodating the
    // divergence components; 

# ifdef REQUIRE
    Require("y and x have same mesh",mesh == x->mesh) ;
    InMethod("FlatField::SetToDivergence(x,vector)") ;
# endif

    int i, j, ncx, dim, icx ;
    FlatField* temp;
  
    ncx = x->nbComponents ;
    dim = mesh->GetDimension() ;  
  
    temp = x->GetWork(1); 

    if (vector){ //vecteur 
	for (i=1 ; i<=dim ; i++){
	    if (i==1)
		SetToGradient(x,i,1,i) ;
	    else{
		temp->SetToGradient(x,i,1,i) ;
		Add(temp,1,1);
	    }  
	}  
    }
    else{// tenseur
	// non-symmetric tensor  
	if (ncx == dim*dim){
	    // Assumes that the tensor components are stored in the 
	    // following order: x11 x12 x13 x21 x22 x23 x31 x32 x33   
	    icx = 1;
	    SetValues(ZERO);
	    for (j=1 ; j<=dim ; j++){ //composante j du vecteur receveur.
		for (i=1 ; i<=dim ; i++){
		    temp->SetToGradient(x,icx++,1,i) ;
		    Add(temp,1,j) ;
		}
	    }
	}
	// symmetric tensor 
	else{
	    if (dim == 2){
		// Assumes that the tensor components are stored in the 
		// following order: x11 x21 x22
		SetToGradient(x,1,1,1) ;
		temp->SetToGradient(x,2,1,2) ;
		Add(temp,1,1) ;
		SetToGradient(x,2,2,1) ;
		temp->SetToGradient(x,3,1,2) ;
		Add(temp,1,2) ;
	    }
	    else{
		// Assumes that the tensor components are stored in the 
		// following order: x11 x21 x31 x22 x32 x33    
		SetToGradient(x,1,1,1) ;
		temp->SetToGradient(x,2,1,2) ;
		Add(temp,1,1) ;
		temp->SetToGradient(x,3,1,3) ;
		Add(temp,1,1) ;
		SetToGradient(x,2,2,1) ;
		temp->SetToGradient(x,4,1,2) ;
		Add(temp,1,2) ;
		temp->SetToGradient(x,5,1,3) ;
		Add(temp,1,2) ;
		SetToGradient(x,3,3,1) ;
		temp->SetToGradient(x,5,1,2) ;
		Add(temp,1,3) ;
		temp->SetToGradient(x,6,1,3) ;
		Add(temp,1,3) ;
	    }
	} 
    }
    x->Retrieve(temp) ;    
}

void FlatField :: SetToDyadicProduct (FlatField* x)
{
    // Initializes the values of the receiver to the dyadic product of 'x x',
    // x is supposed to be a vector.
  
#ifdef REQUIRE
    Require("Valid number of components", nbComponents >= 2 && nbComponents <= 3) ;
    InMethod("FlatField::SetToDyadicProduct(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX ;
    int             k, nc ;

    nc = x->GetNbComponents();

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(x->index) ;
        if (nc == 2){
            eY->CopyFrom(eX,1,1);
            eY->Multiply(eX,1,1);
            eY->CopyFrom(eX,1,2);
            eY->Multiply(eX,2,2);
            eY->CopyFrom(eX,2,3);
            eY->Multiply(eX,2,3); 
        }
        else { // nc == 3
            eY->CopyFrom(eX,1,1);
            eY->Multiply(eX,1,1);
            eY->CopyFrom(eX,1,2);
            eY->Multiply(eX,2,2);
            eY->CopyFrom(eX,1,3);
            eY->Multiply(eX,3,3); 
            eY->CopyFrom(eX,2,4);
            eY->Multiply(eX,2,4);
            eY->CopyFrom(eX,2,5);
            eY->Multiply(eX,3,5);
            eY->CopyFrom(eX,3,6);
            eY->Multiply(eX,3,6);
        } 
    }

}

void FlatField :: SetToGradient (FlatField* x)
{
    // y = dx/dX.
    // Initializes the receiver 'y' to the gradient (with respect to the 
    // coordinates (X,Y,Z)) of 'x'.
    // All components are concerned.
    // If the parent elements of 'y' differ from those of 'x', re-interpolation
    // from 'x' is provided.
    // The receiver is expected to have enough components for accomodating the
    // gradient components; these are stored according to the following 
    // convention: all partial derivatives of the first component are stored 
    // first.

# ifdef REQUIRE
    Require("y and x have same mesh",mesh == x->mesh) ;
    Require("y has correct number of components", nbComponents == x->nbComponents * mesh->GetDimension()) ;
    InMethod("FlatField::SetToGradient(x)") ;
# endif

    int i, j, icy, ncx, dim ;

    icy = 1 ;
    ncx = x->nbComponents ;
    dim = mesh->GetDimension() ;
    for (j=1 ; j<=ncx ; j++)
	for (i=1 ; i<=dim ; i++)
	    SetToGradient(x,j,icy++,i) ;
}


void FlatField :: SetToGradient (FlatField* x, int icx, int icy, int dir)
{
    // y = dx/dX.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th coordinate (X, Y or Z)) of the 'icx'-th
    // component of 'x'.
    // For ex : in 2D : du/dX = du/dr dr/dX + du/ds ds/dX  (where u stands for x)
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', interpolation from 'x' to 'y' is provided.

# ifdef REQUIRE
    Require("y and x have same mesh",mesh == x->mesh) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("valid argument 'dir'", dir >= 1 && dir <= mesh->GetDimension()) ;
    InMethod("FlatField::SetToGradient(x,icx,icy,dir)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eInvJacobian ;
    int             k, indexInvJacobian ;

    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eY->SetToGradient(eX,icx,icy,dir,eInvJacobian) ;
    }
}

void FlatField :: SetToGradientT (FlatField* x, int icx, int icy, int dir)
{
    // y = dx/dX.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th coordinate (X, Y or Z)) of the 'icx'-th
    // component of 'x'.
    // For ex : in 2D : du/dX = du/dr dr/dX + du/ds ds/dX.
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', transposed interpolation from 'x' to 'y' is provided.

# ifdef REQUIRE
    Require("y and x have same mesh",mesh == x->mesh) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("valid argument 'dir'", dir >= 1 && dir <= mesh->GetDimension()) ;
    InMethod("FlatField::SetToGradientT(x,icx,icy,dir)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eInvJacobian ;
    int             k, indexInvJacobian ;

    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eY->SetToGradientT(eX,icx,icy,dir,eInvJacobian) ;
    }
}


void FlatField :: SetToInvJacobian ()
{
    // Initializes the values of the receiver to the partial derivatives of the
    // local coordinates (r,s,t) of every element with respect to the cartesian
    // coordinates (X,Y,Z).
    //
    // The order in which the components are stored is the following:
    //   . 1D: (dr/dX);
    //   . 2D: (dr/dX, dr/dY, ds/dX, ds/dY);
    //   . 3D: (dr/dX, dr/dY, dr/dZ, ds/dX, ds/dY, ds/dZ, dt/dX, dt/dY, dt/dZ).
    //   
    // This method is intended to be used by the inverse jacobian matrix of the
    // meshes, for initializing its values.
    //
    // If the FAST_GRADIENT compilation option is defined, those of the
    // components of the receiver that are close to 0 are set to 0 exactly. The 
    // effect is to speed up subsequent computations of gradients for straight
    // elements (typically, rectangular or parallelipipedic).
    //
    // Note: this method is implemented only for a receiver defined a mesh for
    //    which dim = number of coordinates. For example, the mesh cannot be the
    //    1D skeleton of a 2D mesh (dim = 1, number of coordinates = 2).

# ifdef REQUIRE
    Require("valid mesh type", mesh->GetDimension() == mesh->GetNbCoordinates());
    Require("valid number of components", 
	    nbComponents == mesh->GetDimension() * mesh->GetDimension()) ;
    InMethod("FlatField::SetToInvJacobian()") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eCoord, *eJacobian ;
    int             k, indexCoord, indexJacobian ;

    indexCoord    = mesh->GetCoordinates()->GetIndex() ;
    indexJacobian = mesh->GetJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY        = elem->GetElementaryField(index) ;
        eCoord    = elem->GetElementaryField(indexCoord) ;
        eJacobian = elem->GetElementaryField(indexJacobian) ;
        eY->SetToInvJacobian(eCoord,eJacobian) ;
    }
}


void FlatField :: SetToJacobian ()
{
    // Initializes the values of the receiver to the jacobian, i.e., to the
    // absolute value of the determinant of the matrix of the partial derivatives
    // of the cartesian coordinates (X,Y,Z) with respect to the local coordinates
    // (r,s,t) of every element.
    //
    // This method is intended to be used by the jacobian field of the meshes,
    // for initializing its values.

# ifdef REQUIRE
    Require("the field has 1 component", nbComponents == 1) ;
    InMethod("FlatField::SetToJacobian()") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eCoord ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY     = elem->GetElementaryField(index) ;
        eCoord = elem->GetElementaryField(mesh->GetCoordinates()->GetIndex()) ;
        eY->SetToJacobian(eCoord) ;
    }
}


void FlatField :: SetToLaplacian (FlatField* x)
{
    // Initializes the receiver 'y' to the (strong) Laplacian of 'x'.
    // If the parent elements of 'y' differ from those of 'x', re-interpolation
    // from 'x' is provided.
    //
    // So far, the only use of this method is for evaluating a STRONG Laplacian
    // operator, as required in the evaluation of the k_i coefficients of a 
    // Runge-Kutta time integration scheme.

# ifdef REQUIRE
    Require("y and x have same mesh",mesh == x->mesh) ;
    Require("y and x have same nb components", 
	    nbComponents == x->GetNbComponents()) ;
    InMethod("FlatField::SetToLaplacian(x)") ;
# endif

    FlatField *grad, *work ;
    int       i, j, icgrad, dim ;

    // 1. get gradient

    dim  = mesh->GetDimension() ;
    grad = x->GetWork(nbComponents*dim) ;
    grad->SetToGradient(x) ;

    // 2. get divergence of gradient

    SetValues(ZERO) ;
    work   = GetWork(1) ;
    icgrad = 1 ;
    for (j=1 ; j<=nbComponents ; j++)
	for (i=1 ; i<=dim ; i++) {
	    work->SetToGradient(grad,icgrad++,1,i) ;
	    Add(work,1,j) ;
	}

    x->Retrieve(grad) ;
    Retrieve(work) ;
}


void FlatField :: SetToNormal ()
{
    // Initializes the values of the receiver to the normal vector "n", which is
    // a function of the partial derivatives of the cartesian coordinates (X,Y,Z)
    // with respect to the local coordinates (r,s,t) of the element.
    // At every dof, n has euclidian norm equal to 1.
    //
    // This method is typically used for initializing the values of the normal 
    // vector of a skin mesh involved in a Neumann condition.
    //
    // See more information in method ElementaryField::SetToNormal.

# ifdef REQUIRE
    Require("valid nb of components", nbComponents == mesh->GetNbCoordinates()) ;
    Require("valid type of mesh",
	    (mesh->GetDimension() == 1 && mesh->GetNbCoordinates() == 2)
	    || (mesh->GetDimension() == 2 && mesh->GetNbCoordinates() == 3)) ;
    InMethod("FlatField::SetToNormal()") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eCoord ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY     = elem->GetElementaryField(index) ;
        eCoord = elem->GetElementaryField(mesh->GetCoordinates()->GetIndex()) ;
        ElementaryField *jac = eY->SetToNormal(eCoord) ;
        delete jac;  // Required to avoid memory leak.
    }
}


void FlatField :: SetToNormByElement (FlatField* x, string type)
{
    // The receiver is assumed to possess 1 value per elementary field, i.e., to
    // be interpolated by the GL 0 rule everywhere.
    // Initializes, for every component of for every element field, that unique
    // value to the norm 'type' (e.g., "L2", "H1") of 'x) of all components of
    // the corresponding elementary field of 'x'.

# ifdef REQUIRE
    Require("same mesh", x->mesh == mesh) ;
    Require("one-component field", nbComponents == 1) ;
    Require("valid argument 'type'", !strcmp(type,"L2") ||
	    !strcmp(type,"H1") ||
	    !strcmp(type,"euclidian") ||
	    !strcmp(type,"infinite")) ;
    InMethod("FlatField::SetToNormByElement(x,type)") ;
# endif

    FlatField *temp ;
    Element   *elem, *elem2 ;
    real      norm ;
    int       i, k, kk, n ;

    temp = x->GetWork() ;
  
    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;

	// set values of temp to 0 everywhere, except on 'elem1' (set to receiver)
	temp->CopyFrom(x) ;
	for (kk=1 ; kk<=n ; kk++) {
	    elem2 = temp->mesh->At(kk) ;
	    if (elem2->GetProcessorNumber() == Mpi_rank && elem2 != elem)
		for (i=1 ; i<=nbComponents ; i++)
		    elem2->GetElementaryField(temp->index)->SetValues(ZERO,i) ;
	}

	// set norm of (receiver on 'elem1') to norm of (temp)
	if (! strcmp(type.c_str(),"L2"))
	    norm = temp->GetL2Norm() ;
	else if (! strcmp(type.c_str(),"H1"))
	    norm = temp->GetH1Norm() ;
		 else if (! strcmp(type.c_str(),"euclidian"))
	    norm = temp->GetEuclidianNorm() ;
		 else if (! strcmp(type.c_str(),"infinite"))
	    norm = temp->GetInfiniteNorm() ;
	if (elem->GetProcessorNumber() == Mpi_rank)
	    elem->GetElementaryField(index)->SetValues(norm,1) ;           // 1 value 
    }

    x->Retrieve(temp) ;
}

void FlatField :: SetToRotational (FlatField* x)
{
    // y = grad ^ x.
    // Initializes the receiver 'y' to the (strong) rotational of x.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("correct nb components", 
	    (mesh->GetDimension()==3 && x->nbComponents==3 && nbComponents==3)
	    || (mesh->GetDimension()==2 && x->nbComponents==2 && nbComponents==2)
	    || (mesh->GetDimension()==2 && x->nbComponents==2 && nbComponents==1)); // bcos of stream function ...
    InMethod("FlatField::SetToRotational(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eY2, *eInvJac ;
    int             k, indexInvJacobian  ;

    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY      = elem->GetElementaryField(index) ;
        eX      = elem->GetElementaryField(x->index) ;
        eInvJac = elem->GetElementaryField(indexInvJacobian) ;
        eY2     = eY->GetWork() ;

        if (mesh->GetDimension() == 3) {
            // y1
            eY ->SetToGradient(eX,3,1,2,eInvJac) ;     //  dx3 / dX2
            eY2->SetToGradient(eX,2,1,3,eInvJac) ;     //  dx2 / dX3
            eY ->Subtract(eY2,1,1) ;

            // y2
            eY ->SetToGradient(eX,1,2,3,eInvJac) ;     //  dx1 / dX3
            eY2->SetToGradient(eX,3,1,1,eInvJac) ;     //  dx3 / dX1
            eY ->Subtract(eY2,1,2) ;

            // y3
            eY ->SetToGradient(eX,2,3,1,eInvJac) ;     //  dx2 / dX1
            eY2->SetToGradient(eX,1,1,2,eInvJac) ;     //  dx1 / dX2
            eY ->Subtract(eY2,1,3) ;
        }

        else {        // dim = 2
            eY ->SetToGradient(eX,2,1,1,eInvJac) ;     //  dx2 / dX1
            eY2->SetToGradient(eX,1,1,2,eInvJac) ;     //  dx1 / dX2
            eY ->Subtract(eY2,1,1) ;
        }

        eY->Retrieve(eY2) ;
    }
}

void FlatField :: SetToRotRot (FlatField* x)
{
    // y = grad ^ grad ^ x.= (grad div - Laplacian) x
    // Initializes the receiver 'y' to the (strong) rotational of rotational of  x.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and x are not the same object", this != x) ;
    InMethod("FlatField::SetToRotRot(x)") ;
# endif
  
    FlatField *workx, *worky;

    workx = x->GetWork(1) ;
    worky = GetWork() ;

    workx->SetToDivergence (x, true);
    SetToGradient(workx);
    worky->SetToLaplacian(x);

    Subtract(worky);

    x->Retrieve(workx);
    Retrieve(worky);
}

void FlatField :: SetToRotRot_old (FlatField* x)
{
    // y = grad ^ grad ^ x.
    // Initializes the receiver 'y' to the (strong) rotational of rotational of  x.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("correct nb components", 
	    (mesh->GetDimension()==3 && x->nbComponents==3 && nbComponents==3)
	    || (mesh->GetDimension()==2 && x->nbComponents==2 && nbComponents==2));
    InMethod("FlatField::SetToRotRot(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *rotY, *eY2,*eInvJac ;
    int             k, indexInvJacobian  ;

    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY      = elem->GetElementaryField(index) ;
        eX      = elem->GetElementaryField(x->index) ;
        eInvJac = elem->GetElementaryField(indexInvJacobian) ;
        rotY    = eY->GetWork() ;
        eY2     = eY->GetWork() ;

        if (mesh->GetDimension() == 3) {
            // compute rotational
            // y1
            rotY ->SetToGradient(eX,3,1,2,eInvJac) ;     //  dx3 / dX2
            eY2->SetToGradient(eX,2,1,3,eInvJac) ;     //  dx2 / dX3
            rotY ->Subtract(eY2,1,1) ;//  dx3 / dX2 - dx2 / dX3

            // y2
            rotY ->SetToGradient(eX,1,2,3,eInvJac) ;     //  dx1 / dX3
            eY2->SetToGradient(eX,3,1,1,eInvJac) ;     //  dx3 / dX1
            rotY ->Subtract(eY2,1,2) ; // dx1 / dX3 - dx3 / dX1

            // y3
            rotY ->SetToGradient(eX,2,3,1,eInvJac) ;     //  dx2 / dX1
            eY2->SetToGradient(eX,1,1,2,eInvJac) ;     //  dx1 / dX2
            rotY ->Subtract(eY2,1,3) ;// dx2 / dX1 - dx1 / dX2

            // compute rotational of rotational
            // rot rot 1
            eY->SetToGradient(rotY,3,1,2,eInvJac) ; 
            eY2->SetToGradient(rotY,2,1,3,eInvJac) ;
            eY->Subtract(eY2,1,1);

            // rot rot 2
            eY->SetToGradient(rotY,1,2,3,eInvJac) ; 
            eY2->SetToGradient(rotY,3,1,1,eInvJac) ;
            eY->Subtract(eY2,1,2);

            // rot rot 3
            eY->SetToGradient(rotY,2,3,1,eInvJac) ; 
            eY2->SetToGradient(rotY,1,1,2,eInvJac) ;
            eY->Subtract(eY2,1,3);

        }

        else {        // dim = 2
            rotY ->SetToGradient(eX,2,1,1,eInvJac) ;     //  dx2 / dX1
            eY2 ->SetToGradient(eX,1,1,2,eInvJac) ;     //  dx1 / dX2
            rotY ->Subtract(eY2,1,1) ;// dx2 / dX1 - dx1 / dX2

            eY->SetToGradient(rotY,1,1,2,eInvJac) ; 
            rotY->Multiply(ONE_MINUS,1);
            eY->SetToGradient(rotY,1,2,1,eInvJac) ; 
        }

        eY->Retrieve(rotY) ;
        eY->Retrieve(eY2) ;
    }
}


void FlatField :: SetToWeakConvection (FlatField* x1, FlatField* x2, 
                                       FlatField* rule)
{
    // y = (x1 grad x2, psi).
    // Initializes the receiver 'y' to the weak convection of 'x1' and 'x2',
    // using 'rule' for performing numerical integrations.
    //
    // The convective form (as opposed to the divergence form or to the skew-
    // symmetric form) is used.
    //
    // At method's exit, the values of 'rule' contain garbage.

# ifdef REQUIRE
    Require("y and x1 have same mesh", mesh == x1->mesh) ;
    Require("y and x2 have same mesh", mesh == x2->mesh) ;
    Require("y and rule have same mesh", mesh == rule->mesh) ;
    Require("y and x1 are not the same object", this != x1) ;
    Require("y and x2 are not the same object", this != x2) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("x1 and rule are not the same object", x1 != rule) ;
    Require("x2 and rule are not the same object", x2 != rule) ;
    Require("x1 and y have same nb components",x1->nbComponents == nbComponents);
    Require("x2 and y have same nb components",x2->nbComponents == nbComponents);
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakConvection(x1,x2,rule)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eY2, *eX1, *eX2, *eRule, *eRule2, *eJacobian,
	*eInvJacobian ;
    int             i, k, icomp, indexJacobian, indexInvJacobian ;

    indexJacobian    = mesh->GetJacobian()->GetIndex() ;
    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;
 
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX1          = elem->GetElementaryField(x1->index) ;
        eX2          = elem->GetElementaryField(x2->index) ;
        eRule        = elem->GetElementaryField(rule->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eY2          = eY   ->GetWork() ;
        eJacobian    = eRule->GetWork() ;
        eRule2       = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian),
                                       1,1) ;

        eY->SetValues(ZERO) ;
        for (icomp=1 ; icomp<=nbComponents ; icomp++)
            for (i=1 ; i<=mesh->GetDimension() ; i++) {
                eRule ->SetToGradient(eX2,icomp,1,i,eInvJacobian) ;
                eRule2->CopyInterpolateFrom(eX1,i,1) ;
                eRule ->Multiply(eRule2) ;
                eRule ->Multiply(eJacobian) ;
                eRule ->MultiplyByWeights(1) ;
                eY2   ->CopyInterpolateTFrom(eRule,1,1) ;
                eY    ->Add(eY2,1,icomp) ;
            }
        eY->Retrieve(eY2) ;
        eRule->Retrieve(eJacobian) ;
        eRule->Retrieve(eRule2) ;
    }
}

void FlatField :: SetToWeakDivergence (FlatField* x, FlatField* rule)
{
    // y = (div x, psi).
    // Initializes the receiver 'y' to the weak divergence of 'x', using 'rule'
    // for performing numerical integrations.
    //
    // At method's exit, the values of 'rule' contain garbage.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and rule have same mesh", mesh == rule->mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("x and rule are not the same object", x != rule) ;
    Require("nb components of x = dimension of mesh", 
	    x->nbComponents == mesh->GetDimension()) ;
    Require("y has 1 component", nbComponents == 1) ;
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakDivergence(x,rule)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eWorkY, *eRule, *eJacobian, *eInvJacobian ;
    int             i, k, indexJacobian, indexInvJacobian ;

    indexJacobian    = mesh->GetJacobian()->GetIndex() ;
    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eRule        = elem->GetElementaryField(rule->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eWorkY       = eY->GetWork() ;
        eJacobian    = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian), 1,1) ;

        eY->SetValues(ZERO) ;
        for (i=1 ; i<=x->nbComponents ; i++) {
            eRule ->SetToGradient(eX,i,1,i,eInvJacobian) ;
            eRule ->Multiply(eJacobian) ;
            eRule ->MultiplyByWeights(1) ;
            eWorkY->CopyInterpolateTFrom(eRule,1,1) ;
            eY    ->Add(eWorkY,1,1) ;
        }

        eY->Retrieve(eWorkY) ;
        eRule->Retrieve(eJacobian) ;
    }
}

void FlatField :: SetToWeakGradient (FlatField* x, FlatField* rule)
{
    // y = - (x, grad psi), i.e., y = - Dt x.
    // Initializes the receiver 'y' to the weak gradient of 'x', using 'rule' for
    // performing numerical integrations.
    // 
    // The components of 'y' are ordered according to the following convention:
    // let 'dim' stand for the mesh's dimension (= number of coordinates); then
    // the 'dim' components associated with the first component of 'x' are stored
    // first; the 'dim' components associated with the second component of 'x'
    // are stored next, etc.
    //
    // At method's exit, the values of 'rule' contain garbage.
    //
    // Note: this operator is equal to "minus" the transposition of the weak
    //       divergence operator.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and rule have same mesh", rule->mesh == mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("x and rule are not the same object", x != rule) ;
    Require("nb components of y = nb components of x * dimension of mesh",
	    nbComponents == x->nbComponents * mesh->GetDimension()) ;
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakGradient(x,rule)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eRule, *eJacobian, *eInvJacobian ;
    int             k, icy, icomp, dir, indexJacobian, indexInvJacobian ;

    indexJacobian    = mesh->GetJacobian()->GetIndex() ;
    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eRule        = elem->GetElementaryField(rule->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eJacobian    = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian), 1,1) ;
        icy = 0 ;
        for (icomp=1 ; icomp<=x->nbComponents ; icomp++)
            for (dir=1 ; dir<=mesh->GetDimension() ; dir++) {
                icy++ ;
                eRule->CopyInterpolateFrom(eX,icomp,1) ;
                eRule->MultiplyAndSwitchSigns(eJacobian,1,1) ;
                eRule->MultiplyByWeights(1) ;
                eY   ->SetToGradientT(eRule,1,icy,dir,eInvJacobian) ;
            }
        eRule->Retrieve(eJacobian) ;
    }
}


void FlatField :: SetToWeakIdentity (FlatField* x, FlatField* rule)
{
    // y = (x, psi).
    // Initializes the receiver 'y' to the weak form of 'x', using 'rule' for 
    // performing numerical integrations.
    //
    // At method's exit, the values of 'rule' contain garbage.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and rule have same mesh", rule->mesh == mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("x and rule are not the same object", x != rule) ;
    Require("y and x have same nb components", x->nbComponents == nbComponents) ;
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakIdentity(x,rule)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eRule, *eJacobian ;
    int             i, k, indexJacobian ;

    indexJacobian = mesh->GetJacobian()->GetIndex() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY        = elem->GetElementaryField(index) ;
        eX        = elem->GetElementaryField(x->index) ;
        eRule     = elem->GetElementaryField(rule->index) ;
        eJacobian = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian), 1,1) ;
        for (i=1 ; i<=nbComponents ; i++) {
            eRule->CopyInterpolateFrom(eX,i,1) ;
            eRule->Multiply(eJacobian) ;
            eRule->MultiplyByWeights(1) ;
            eY   ->CopyInterpolateTFrom(eRule,1,i) ;
        }
        eRule->Retrieve(eJacobian) ;
    }
}


void FlatField :: SetToWeakLaplacian (FlatField* x, FlatField* rule)
{
    // y = - (grad x, grad psi).
    // Initializes the receiver 'y' to the weak Laplacian of 'x', using 'rule'
    // for performing numerical integrations.
    //
    // The components of 'y' are ordered according to the following convention:
    // let 'dim' stand for the mesh's dimension (= number of coordinates); then
    // the 'dim' components associated with the first component of 'x' are stored
    // first; the 'dim' components associated with the second component of 'x'
    // are stored next, etc.
    //
    // At method's exit, the values of 'rule' contain garbage.

# ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and rule have same mesh", rule->mesh == mesh) ;
    Require("y and x are not the same object", this != x) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("x and rule are not the same object", x != rule) ;
    Require("y and x have same nb components", nbComponents == x->nbComponents) ;
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakLaplacian(x,rule)") ;
# endif

    Element         *elem ;
    ElementaryField *eY, *eX, *eWorkY, *eRule, *eJacobian, *eInvJacobian ;
    int             i, k, icomp, dim, indexJacobian, indexInvJacobian ;

    indexJacobian    = mesh->GetJacobian()->GetIndex() ;
    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;
    dim              = mesh->GetDimension() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eRule        = elem->GetElementaryField(rule->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eWorkY       = eY->GetWork() ;
        eJacobian    = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian),
                                       1,1) ;
        eY->SetValues(ZERO) ;
        for (icomp=1 ; icomp<=nbComponents ; icomp++)
            for (i=1 ; i<=dim ; i++) {
                eRule ->SetToGradient(eX,icomp,1,i,eInvJacobian) ;
                eRule ->Multiply(eJacobian) ;
                eRule ->MultiplyByWeights(1) ;
                eWorkY->SetToGradientT(eRule,1,1,i,eInvJacobian) ;
                eY    ->Add(eWorkY,ONE_MINUS,1,icomp) ;
            }
        eRule->Retrieve(eJacobian) ;
        eY->Retrieve(eWorkY) ;
    }
}

void FlatField :: SetToZeroMean (int icy)
{
    // Subtracts from the 'icy'-th component the mean value of that component.
    // The mean is computed in the functional sense, i.e., by computing the
    // integral.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::SetToZeroMean(icy)") ;
# endif

    real mean ;

    mean = GetMeanValue(icy) ;

    Add(-mean,icy) ;
}

void FlatField :: SetToWeakTensorDivergence (FlatField* x, FlatField* rule, boolean b)
{
    // y = -(x, grad psi).
    // Initializes the receiver 'y' to the weak divergence of 'x', using 'rule'
    // for performing numerical integrations.
    // x is supposed to be a SYMMETRIC tensor.
    //
    // b must be set to true if integration by part is performed
    //                  false otherwise.
    //
    // At method's exit, the values of 'rule' contain garbage.
    //
    // Used as result for integration by part of the additional viscoleastic term 
    // in the momentum equation.
 
#ifdef REQUIRE
    Require("y and x have same mesh", mesh == x->mesh) ;
    Require("y and rule have same mesh", mesh == rule->mesh) ;
    Require("y and rule are not the same object", this != rule) ;
    Require("y has number of components of mesh dimension", nbComponents == mesh->GetDimension()) ;
    Require("rule has 1 component", rule->nbComponents == 1) ;
    InMethod("FlatField::SetToWeakTensorDivergence(x,rule,boolean)") ;
# endif
 
    Element *elem ;
    ElementaryField *eY, *eX, *eWorkY, *eRule, *eJacobian, *eInvJacobian ;
    int i, j, k, dim, indexJacobian, indexInvJacobian ;
    int iconv;

    indexJacobian    = mesh->GetJacobian()->GetIndex() ;
    indexInvJacobian = mesh->GetInvJacobian()->GetIndex() ;
    dim = mesh->GetDimension() ;
    
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY           = elem->GetElementaryField(index) ;
        eX           = elem->GetElementaryField(x->index) ;
        eRule        = elem->GetElementaryField(rule->index) ;
        eInvJacobian = elem->GetElementaryField(indexInvJacobian) ;
        eWorkY       = eY->GetWork() ;
        eJacobian    = eRule->GetWork() ;
        eJacobian->CopyInterpolateFrom(elem->GetElementaryField(indexJacobian),1,1) ;

        eY->SetValues(ZERO) ;
        for (j=1 ; j<=dim ; j++) {
            for (i=1 ; i<=dim ; i++) {
    
                //        iconv = (int) iD->At(j,i) ;
                iconv = ( j - 1 ) * dim + i; // j=colonne, i=lignes.

                if (b){
                    // Terms to be used if integration by parts in the momentum equation
                    eRule->CopyInterpolateFrom(eX,iconv,1) ;
                    eRule->MultiplyAndSwitchSigns(eJacobian,1,1) ;
                    eRule->MultiplyByWeights(1) ;
                    eWorkY->SetToGradientT(eRule,1,1,i,eInvJacobian) ;
                    eY->Add(eWorkY,1,j) ;
                }
                else{  
                    // Terms to be used if no integration by parts in the momentum equation
                    eRule->SetToGradient(eX,iconv,1,i,eInvJacobian) ;
                    eRule->Multiply(eJacobian) ;
                    eRule->MultiplyByWeights(1) ;
                    eWorkY->CopyInterpolateTFrom(eRule,1,1) ;
                    eY->Add(eWorkY,1,j) ;
                }
            }
        }
        eY->Retrieve(eWorkY) ;
        eRule->Retrieve(eJacobian) ;
    }
}  

void FlatField :: SetToZeroAlgebraicMean (int icy)
{
    // Subtracts from the 'icy'-th component the mean value of that component.
    // The mean is computed in the algebraic sense, i.e., by summing (without
    // weighting) the values.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >=1 && icy <= nbComponents) ;
    InMethod("FlatField::SetToZeroAlgebraicMean(icy)") ;
# endif

    Element    *elem ;
    RealVector *comp ;
    real       sum, mean ; 
    int        k, m, size ;

    sum = ZERO ;
    m   = 0 ;
    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        comp = elem->GetElementaryField(index)->GetComponent(icy) ;
        sum += comp->SumValues() ;
        m   += comp->GetSize() ;
    }

    mean = Mpi_Sum(sum) ;
    size = Mpi_Sum(m) ;
    mean /= size ;

    Add(-mean,icy) ;
}
  
void FlatField :: SetValues (real val)
{
    // Sets every value of the receiver to 'val'.

    Element         *elem ;
    ElementaryField *eY ;
    int             k, ic ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        for (ic=1 ; ic<=nbComponents ; ic++)
            eY->SetValues(val,ic) ;
    }
}


void FlatField :: SetValues (real val, int icy)
{
    // Sets every value of the 'icy'-th component of the receiver to 'val'.

    Element *elem ;
    int     k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        elem->GetElementaryField(index)->SetValues(val,icy) ;
    }
}

void FlatField :: ShiftPolynomialDegree (int shift)
{
    // Increments by 'shift' the polynomial degree of every elementary field of
    // the receiver.
    // The current values are lost and reinitialized to 0.

    Element *elem ;
    int     k, n ;

    n = mesh->GetSize() ;
    for (k=1 ; k<=n ; k++) {
	elem = mesh->At(k) ;
	elem->GetElementaryField(index)->ShiftPolynomialDegree(shift) ;
    }

    nbDofs         = UNINITIALIZED_INTEGER ;
    nbInteriorDofs = UNINITIALIZED_INTEGER ;
}

real FlatField :: SumValues ()
{
    // Returns the sum of all values of all components of the receiver.

    Element *elem ;
    real    sum, answer ;
    int     k;

    sum = ZERO ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        sum += elem->GetElementaryField(index)->SumValues() ;
    }

    answer = Mpi_Sum(sum) ;

    return answer ;
}


void FlatField :: Square()
{  
    Element         *elem ;
    ElementaryField *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eX = elem->GetElementaryField(index) ;
        eX->Square();
    }
}      

void FlatField :: Sqrt()
{ 

    Element         *elem ;
    ElementaryField *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eX = elem->GetElementaryField(index) ;
        eX->Sqrt();
    }
}


void FlatField :: Power(real x)
{  
    //MAH 04.2006

    Element         *elem ;
    ElementaryField *eX ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eX = elem->GetElementaryField(index) ;
        eX->Power(x);
    }
}


void FlatField :: ScalarProductWithTensors (FlatField* x1, FlatField* x2)
{   // MAH 12.2005
    // initialise le receveur avec le produit tensoriel de x1 x2
    
#ifdef REQUIRE
    Require("x1 is a flat field", x1->IsFlatField()) ;
    Require("x2 is a flat field", x2->IsFlatField()) ;
    Require("y and x1 have same mesh", mesh == x1->GetMesh()) ;
    Require("y and x2 have same mesh", mesh == x2->GetMesh()) ;
    Require("y and x2 have same interpolation", HasSameInterpolationAs(x2)) ;
    Require("y and x1 have same interpolation", HasSameInterpolationAs(x1)) ;
    Require("valid argument 'icx'", nbComponents==1) ;
    InMethod("FlatField::(x1,x2)") ;
# endif

    Element         *elem ;
    ElementaryField *eX1,*eX2,*eY ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX1 = elem->GetElementaryField(x1->GetIndex()) ;
        eX2 = elem->GetElementaryField(x2->GetIndex()) ;
        eY->SetValues(ZERO) ;//initialise le receveur avec le produit tensoriel.

        eY->ScalarProductWithTensors(eX1, eX2);
    }
}

void FlatField :: SetToTrace (FlatField* x)
{   // MAH 07.12.2005
    // Receiver gets the trace of x
    // x must be a tensor of order 2 with 9 components
    
#ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y and x have same mesh", mesh == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", nbComponents==1) ;
    Require("",x->GetNbComponents()==9);
    InMethod("FlatField::(x)") ;
# endif

    Element         *elem ;
    ElementaryField *eX,*eY ;
    int             k;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        eY = elem->GetElementaryField(index) ;
        eX = elem->GetElementaryField(x->GetIndex()) ;
        eY->SetValues(ZERO) ;

        eY->SetToTrace(eX);
    }
}


//______________________________ MortaredField ________________________________


MortaredField :: MortaredField ()
    : Field()
{
    // Protected constructor, to be used only internally, e.g. in methods for
    // duplicating mortared fields.

    main   = NULL ;
    mortar = NULL ;
}


MortaredField :: MortaredField (int ncomp, Mesh* m, ParentElement* p)
    : Field ()
{
    // Constructor. Initializes the receiver to a mortared field whose main field
    // is defined on 'm', with 'ncomp' components. Every elementary field of the
    // main field has 'p' as parent element in every direction. 
    // The values are initialized to 0.

# ifdef REQUIRE
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("MortaredField::MortaredField(ncomp,m,p)") ;
# endif

    main = new FlatField(ncomp,m,p) ;
    CreateMortar() ;
}


MortaredField :: MortaredField (int ncomp, Mesh* m, Vector<ParentElement*>* vp)
    : Field ()
{
    // Constructor. Initializes the receiver to a mortared field whose main field
    // is defined on 'm', with 'ncomp' components. For every elementary field, 
    // the parent element vp(i) is used as parent element in the i-th direction.
    // The values are initialized to 0.
    // Remark: every elementary field will point to a copy of 'vp', not to 'vp'
    //         directly, so 'vp' can be deleted safely.

# ifdef REQUIRE
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("MortaredField::MortaredField(ncomp,m,vp)") ;
# endif

    main = new FlatField(ncomp,m,vp) ; 
    CreateMortar() ;
}


MortaredField :: MortaredField (int ncomp, FlatField* x)
    : Field ()
{
    // Constructor. Initializes the receiver to a mortared field whose main field
    // is defined on the same mesh as 'x', with the same interpolation as 'x', 
    // with 'ncomp' components. 
    // The values are initialized to 0.

    main = x->DuplicateEmpty(ncomp) ;
    CreateMortar() ;
}


MortaredField :: MortaredField (int ncomp, FlatField* x, Mesh* m)
    : Field ()
{
    // Constructor. Initializes the receiver to a mortared field whose main is
    // defined on 'm', with 'ncomp' components. 'm' must be included in the mesh
    // of 'x', and can be of a smaller dimension than the mesh of 'x'.
    // For every element, the parent elements are chosen as the corresponding 
    // parent elements of the spectral element of the mesh of 'x'.
    // The values of the receiver are initialized to 0.
    // See FlatField::FlatField(ncomp,x,m) for additional information.

# ifdef REQUIRE
    Require("mesh 'm' is dispatched", m->IsDispatched()) ;
    InMethod("MortaredField::MortaredField(ncomp,x,m)") ;
# endif

    main = new FlatField(ncomp,x,m) ;
    CreateMortar() ;
}


MortaredField :: ~MortaredField ()
{
    // Destructor.

    unref(main) ;
    unref(mortar) ;
}


void MortaredField :: Add (Field* x)
{
    // y = y + x .
    // Adds 'x' to the receiver 'y'.
    // HP compiler is a shit: this method should be inherited (and inlined).

    Add(x,ONE) ;
}


void MortaredField :: Add (Field* x, real alpha) 
{
    // y = y + alpha x .
    // Adds the product 'alpha*x' to the receiver 'y', recursively.

# ifdef REQUIRE
    Require("x is a mortared field", x->IsMortaredField()) ;
    Require("y and x have same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb of components", 
	    GetNbComponents() == x->GetNbComponents());
    InMethod("MortaredField::Add(x,alpha)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->Add(X->main,alpha) ;                       // upper level

    if (mortar)
	mortar->Add(X->mortar,alpha) ;                 // lower level
}


void MortaredField :: Add (Field* x, int icx, int icy)
{
    // y(icy) = y(icy) + x(icx) .
    // Adds the 'icx'-th component of 'x' to the 'icy'-th component of the
    // receiver 'y'.
    // HP compiler is a shit: this method should be inherited (and inlined).

    Add(x,ONE,icx,icy) ;
}


void MortaredField :: Add (Field*, real, int, int)
{
    // Not implemented yet.

    NotImplemented("MortaredField::Add(x,alpha,icx,icy)") ;
}


void MortaredField :: Add (real, int)
{
    // Not implemented yet.

    NotImplemented("MortaredField::Add(val)") ;
}

void MortaredField :: Assemble ()
{
    // Performs the assembly (i.e., direct-stiffness) operation on the receiver.
    //
    // The values are copied from the main field to the mortar, recursively,
    // using the transposed mortar operator Q(T) wherever non-conformity occurs.

    if (mortar) 
    {
	mortar->main->SetValues(ZERO) ;
	main->AddValuesTo(mortar->main,MORTAR_NC) ;
	mortar->Assemble() ;
    }
}


void MortaredField :: Distribute ()
{
    // Initializes the values of the main field of the receiver.
    // This method is the "de-assembling" counter-part of method Assemble.
    //
    // The values are copied form the mortar to the main field, recursively from
    // dimension 0, using the mortar operator Q wherever non-conformity occurs.
    //
    // This operation does not affect the internal values of the main field; it
    // only affects the values on the contour of every element.

    if (mortar) {
	mortar->Distribute() ;
	main->CopyValuesFrom(mortar->main,MORTAR_NC) ;
    }
}

void MortaredField :: CopyFrom (Field* x)
{
    // y = x.
    // Copies 'x' into 'y', recursively.

# ifdef REQUIRE
    Require("x is a MortaredField", x->IsMortaredField()) ;
    Require("y and x have same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb of components", 
	    GetNbComponents() == x->GetNbComponents());
    InMethod("MortaredField::CopyFrom(x)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->CopyFrom(X->main) ;                            // upper level

    if (mortar)
	mortar->CopyFrom(X->mortar) ;                      // lower level
}


void MortaredField :: CopyFrom (Field* x, int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y', recursively.

# ifdef REQUIRE
    Require("x is a MortaredField", x->IsMortaredField()) ;
    Require("y and x have same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= GetNbComponents()) ;
    InMethod("MortaredField::CopyFrom(x,icx,icy)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->CopyFrom(X->main,icx,icy) ;                    // upper level

    if (mortar)
	mortar->CopyFrom(X->mortar,icx,icy) ;              // lower level
}


void MortaredField :: CopyValuesFrom (FlatField* x, int icx, int icy)
{
    // Initializes the 'icy'-th component of the receiver to the 'icx'-th
    // component of 'x', recursively.
    // Uses normal interpolation mode and mortar conformity mode.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >=1 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >=1 && icy <= GetNbComponents()) ;
    InMethod("MortaredField::CopyValuesFrom(x,icx,icy)") ;
# endif 

    main->CopyValuesFrom(x,icx,icy) ;                      // upper level 

    if (mortar)                            
	mortar->CopyValuesFrom(x,icx,icy) ;                  // lower level
}


void MortaredField :: CreateMortar ()
{
    // Initializes the mortar field of the receiver, recursively until dimension
    // 0.

    Mesh *mesh ;

    mesh = main->GetMesh() ;
    if (mesh->GetDimension() > 0)
	mortar = new MortaredField(GetNbComponents(),main,mesh->GetSkeleton()) ;
    else
	mortar = NULL ;
}

real MortaredField :: Dot (Field* x) 
{
    // dot = y(transposed) . x
    // Returns the scalar product of the receiver 'y' and 'x', by accumulating
    // the dot products on the interior dofs recursively.
    // Involves all components of 'y' and 'x'.

# ifdef REQUIRE
    Require("x is a MortaredField", x->IsMortaredField()) ;
    Require("y and x have same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb of components", 
	    GetNbComponents() == x->GetNbComponents()) ;
    InMethod("MortaredField::Dot(x)") ;
# endif

    MortaredField *X ;
    real          answer ;

    X = (MortaredField*)x ;

    answer = main->InteriorDot(X->main) ;                  // upper level

    if (mortar)
	answer += mortar->Dot(X->mortar) ;                   // lower level

    return answer ;
}


MortaredField* MortaredField :: Duplicate ()
{
    // Returns an exact copy of the receiver.

    MortaredField *answer ;

    answer = DuplicateEmpty(ALL) ;
    answer->CopyFrom(this) ;

    return answer ;
}


MortaredField* MortaredField :: DuplicateEmpty (int ncomp)
{
    // Returns a new mortared field, a copy of the receiver, with 'ncomp' 
    // components and all values initialized to 0.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.

    MortaredField *answer ;

    answer       = new MortaredField() ;
    answer->main = main->DuplicateEmpty(ncomp) ;
    if (mortar)
	answer->mortar = (MortaredField*) mortar->DuplicateEmpty(ncomp) ;

    return answer ;
}


void MortaredField :: EnforceMortarConstraints ()
{
    // Initializes the values of the main field of the receiver, so as to satisfy
    // mortar constraints.
    //
    // First, the values of the main field are copied to the mortar field,
    // recursively. This is done by pointwise re-interpolation.
    // Second, the values are copied back to the main field from the mortar,
    // recursively. This is done by mortar reinterpolation, i.e., using the Q
    // operator.
    //
    // This operation can be seen as: assembly (without summation) + de-assembly.
    // It is rarely used. Its main scope is when the initial values of a field
    // for time-marching problem are given analytically; in that case, mortar
    // constraints are not satisfied at non-conformal interfaces. Before 
    // performing any operation on that field, it is advisable to invoke this
    // method.

    if (mortar) {
	mortar->main->CopyValuesFrom(main,POINTWISE_NC) ;
	mortar->EnforceMortarConstraints() ;
    }

    Distribute() ;
}


real MortaredField :: GetEuclidianNorm ()
{
    // Returns the euclidian (= discrete L2) norm of the receiver, by 
    // accumulating the interior norms of the main field and the mortar field.
    // Considers the receiver as an assembled field.
  
    real normSquared, mortarNorm ;

    normSquared = main->InteriorDot(main) ;            // upper level

    if (mortar) {                                      // lower level
	mortarNorm   = mortar->GetEuclidianNorm() ;
	normSquared += mortarNorm * mortarNorm ;
    }

    return sqrt(normSquared) ;
}


real MortaredField :: GetInfiniteNorm ()
{
    // Returns the L-infinite norm of the receiver, by getting the largest
    // coefficient of the interior main field and the mortar field.
    // Considers the receiver as an assembled field.
  
    real answer ;

    answer = main->GetInteriorInfiniteNorm() ;         // upper level

    if (mortar)                                        // lower level
	answer = max(answer,mortar->GetInfiniteNorm()) ;

    return answer ;
}


int MortaredField :: GetNbDofs ()
{
    // Returns the number of degrees of freedom of the receiver.
    // This number is obtained by accumulating the number of dofs of the main 
    // field and of the mortar field.

    int answer ;

    answer = main->GetNbInteriorDofs() ;               // upper level

    if (mortar)                                        // lower level
	answer += mortar->GetNbDofs() ;

    return answer ;
}


void MortaredField :: GetSetDof (int iop, int idof, real& val)
{
    // If 'iop' = 0, gets into 'val' the value of the 'idof'-th dof.
    // If 'iop' = 1, sets to 'val' the value of the 'idof'-th dof.

# ifdef REQUIRE
    Require("valid dof number", idof >= 1 && idof <= GetNbDofs()) ;
    InMethod("MortaredField::GetSetDof(iop,idof,val)") ;
# endif

    int nbDofMain ;

    nbDofMain = main->GetNbInteriorDofs() ;

    if (idof <= nbDofMain) {
	// the dof is in the main field
	if (mortar)                                // general case
	    main->GetSetInteriorDof(iop,idof,val) ;
	else                                       // particular case of a 0D field
	    main->GetSetDof(iop,idof,val) ;
    }

    else
	// the dof is in the mortar field
	mortar->GetSetDof(iop,idof-nbDofMain,val) ;
}


MortaredField* MortaredField :: GetWork (int ncomp)
{
    // Returns a copy of the receiver, with 'ncomp' components.
    // 'ncomp' = ALL (i.e., = 0) means 'ncomp' = nbComponents.
    // The values of the answer contain garbage.

    MortaredField *answer ;

    answer = new MortaredField() ;

    answer->main = main->GetWork(ncomp) ;
    if (mortar)
	answer->mortar = mortar->GetWork(ncomp) ;

    return answer ;
}


boolean MortaredField :: HasSameInterpolationAs (Field* x)
{
    // Returns 'true' if:
    //  - the receiver is defined on the same mesh as 'x', recursively, and
    //  - on every element the receiver has the same parent elements as 'x',
    //    recursively.

# ifdef REQUIRE
    Require("x is a mortared field", x->IsMortaredField()) ;
    InMethod("MortaredField::HasSameInterpolationAs(x)") ;
# endif

    MortaredField * X ;
 
    X = (MortaredField*)x ;

    if (! main->HasSameInterpolationAs(X->main))          // upper level
	return false ;

    if (mortar)                                           // lower level
	if (! mortar->HasSameInterpolationAs(X->mortar))
	    return false ;

    return true ;
}


boolean MortaredField :: HasSameTypeAs (Field* x)
{
    // Returns 'true' if 'x' is a mortared field, else returns 'false'.

    return x->IsMortaredField() ;
}


void MortaredField :: Inverse ()
{
    // y = 1 / y .
    // Substitutes every value of the receiver by its inverse, recursively.

    main->Inverse() ;                                  // upper level

    if (mortar)
	mortar->Inverse() ;                              // lower level
}


boolean MortaredField :: IsMortaredField ()
{
    // Returns 'true', because the receiver is a mortared field.

    return true ;
}


void MortaredField :: Multiply (real alpha) 
{
    // y = alpha * y .
    // Multiplies the receiver by 'alpha', recurssively.

    main->Multiply(alpha) ;                            // upper level

    if (mortar)
	mortar->Multiply(alpha) ;                        // lower level
}


void MortaredField :: Multiply (Field* x) 
{
    // y = y * x .
    // Multiplies every value of the receiver 'y' by the vis-a-vis values in 'x',
    // recursively.
    // All components are affected.

# ifdef REQUIRE
    Require("x is a mortared field", x->IsMortaredField()) ;
    Require("y and x have the same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components",
	    GetNbComponents() == x->GetNbComponents()) ;
    InMethod("MortaredField::Multiply(x)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->Multiply(X->main) ;                          // upper level

    if (mortar)
	mortar->Multiply(X->mortar) ;                    // lower level
}


void MortaredField :: Multiply (Field* x, int icx, int icy) 
{
    // y = y * x .
    // Multiplies every value of the 'icy'-th component of the receiver 'y' by 
    // the 'icx'-th component of 'x', recursively.

# ifdef REQUIRE
    Require("x is a mortared field", x->IsMortaredField()) ;
    Require("y and x have the same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >= 0 && icx <= x->GetNbComponents()) ;
    Require("valid argument 'icy'", icy >= 0 && icy <= GetNbComponents()) ;
    InMethod("MortaredField::Multiply(x,icx,icy)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->Multiply(X->main,icx,icy) ;                   // upper level

    if (mortar)
	mortar->Multiply(X->mortar,icx,icy) ;             // lower level
}


void MortaredField :: MultiplyAndAdd (real alpha, Field* x, real beta) 
{
    // y = alpha*y + beta*x .
    // Multiplies the receiver 'y' by 'alpha' then adds 'beta * x'.

# ifdef REQUIRE
    Require("x is a MortaredField", x->IsMortaredField()) ;
    Require("y and x have same mesh", GetMesh() == x->GetMesh()) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb of components", 
	    GetNbComponents() == x->GetNbComponents());
    InMethod("MortaredField::MultiplyAndAdd(alpha,x,beta)") ;
# endif

    MortaredField *X ;

    X = (MortaredField*)x ;

    main->MultiplyAndAdd(alpha,X->main,beta) ;              // upper level

    if (mortar)
	mortar->MultiplyAndAdd(alpha,X->mortar,beta) ;        // lower level
}


void MortaredField :: Print ()
{
    // Prints the receiver on standard output. All values are printed.
    // HP compiler is a shit: this method should be inherited.

    Print(true) ;
}

void MortaredField :: Print (boolean showValues)
{
    // Prints the receiver on standard output. Prints the values only if
    // 'showValues' is set to 'true'.
    // Only the main field is printed.

  if (!name.empty())
    Printf("Mortared field '%s' whose main field is:\n",name.c_str()) ;
    else
	Printf("Mortared field whose main field is:\n") ;

    main->Print(showValues) ;

    if (mortar)
	Printf("mortar not printed\n\n") ;
}
  

void MortaredField :: PrintComponent (int icy)
{
    // Prints on standard output the 'icy'-th component of the receiver.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 0 && icy <= GetNbComponents()) ;
    InMethod("MortaredField::Print(icy)") ;
# endif
  
    if (!name.empty())
      Printf("Mortared field '%s' whose main field is:\n",name.c_str()) ;
    else
	Printf("Mortared field whose main field is:\n") ;

    main->PrintComponent(icy) ;

    if (mortar) {
	Printf("mortar is: ") ;
	mortar->PrintComponent(icy) ;
    }
}
  

void MortaredField :: PrintWithMortar ()
{
    // Prints the receiver on standard output.
    // The mortar is also printed.

  if (!name.empty())
    Printf("Mortared field '%s' whose main field is:\n",name.c_str()) ;
    else
	Printf("Mortared field whose main field is:\n") ;

    main->Print(true) ;

    if (mortar) {
	Printf("mortar is: ") ;
	mortar->PrintWithMortar() ;
    }
}
  

void MortaredField :: Retrieve (Field* w)
{
    // Makes 'w' available again as a work field identical to the receiver.
    // ('w' may have a different number of components from the receiver).
    // With respect to reference counting, invoking this method is tantamount
    // to calling "unref(w)".

# ifdef REQUIRE
    Require("'w' is a mortared field", w->IsMortaredField()) ;
    Require("same interpolation",  HasSameInterpolationAs(w)) ;
    InMethod("MortaredField::Retrieve(w)") ;
# endif

    MortaredField *ww ;

    ww = (MortaredField*) w ;

    main->Retrieve(ww->main) ;
    if (mortar)
	mortar->Retrieve(ww->mortar) ;

    ww->main   = NULL ;
    ww->mortar = NULL ;
    delete ww ;
}

void MortaredField :: SetParentElements (ParentElement* p) 
{
    // Changes the discretization rule in every elementary field and in every
    // direction to 'p'.
    // Applies also to the mortar.
    // All values are lost and reinitialized to 0.

    main->SetParentElements(p) ; 

    if (mortar)
	mortar->SetParentElements(p) ;
}


void MortaredField :: SetParentElements (ParentElement* p, int dir) 
{
    // Changes the discretization rule in every elementary field to 'p', in the
    // 'dir'-th direction.
    // Applies also to the mortar.
    // All values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("valid argument 'dir'", dir>=1 && dir <= GetMesh()->GetDimension()) ;
    InMethod("MortaredField::SetParentElements(p,dir)") ;
# endif

    main->SetParentElements(p,dir) ; 

    if (mortar)
	mortar->SetParentElements(p,dir) ;
}


void MortaredField :: SetParentElement (ParentElement* p, int dir, 
                                        Element* elem)
{
    // Changes the discretization rule in the elementary field of 'elem', in the
    // 'dir'-th direction, to 'p'.
    // Applies also to the mortar.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("'elem' belongs to the mesh", GetMesh()->Has(elem)) ;
    Require("valid argument 'dir'", dir>=1 && dir<=GetMesh()->GetDimension()) ;
    InMethod("MortaredField::SetParentElements(p,dir,elem)") ;
# endif

    main->SetParentElement(p,dir,elem) ;

    if (mortar)
	mortar->SetParentElement(p,dir,elem) ;
}


void MortaredField :: SetPolynomialDegree (int n)
{
    // Sets to 'n' the polynomial degree of all elementary fields of the 
    // receiver, in all directions.
    // Applies also to the mortar.
    // The type of the collocation rule (GL, GLL, etc) need not be the same in
    // all elementary fields, nor in all directions within every elementary 
    // field.
    // The current values are lost and reinitialized to 0.

    main->SetPolynomialDegree(n) ;

    delete mortar ;
    CreateMortar() ;
}


void MortaredField :: SetPolynomialDegree (int n, int dir)
{
    // Sets to 'n' the polynomial degree of all elementary fields of the 
    // receiver, in the 'dir'-th directions. 
    // Applies also to the mortar.
    // The type of the collocation rule (GL, GLL, etc) need not be the same in
    // all elementary fields, nor in all directions within every elementary 
    // field.
    // The current values are lost and reinitialized to 0.

    main->SetPolynomialDegree(n,dir) ;

    delete mortar ;
    CreateMortar() ;
}


void MortaredField :: SetPolynomialDegree (int n, int dir, Element* elem)
{
    // Sets to 'n' the polynomial degree of the elementary field of 'elem', in
    // the 'dir'-th directions. 
    // Applies also to the mortar.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("'elem' belongs to the mesh", GetMesh()->Has(elem)) ;
    InMethod("MortaredField::SetPolynomialDegree(n,dir,elem)") ;
# endif

    main->SetPolynomialDegree(n,dir,elem) ;

    delete mortar ;
    CreateMortar() ;
}


void MortaredField :: SetTo (Function* f, real t)
{
    // Initializes the values of the receiver to the values of 'f'. Any dof with
    // coordinates (x,y,z) is initialized to f(x,y,z,t).
    // All components are initialized and thus become identical.

    main->SetTo(f,t) ;

    if (mortar)
	mortar->SetTo(f,t) ;
}


void MortaredField :: SetTo (Function* f, int ic, int icf, real t)
{
    // Initializes the values of the 'ic'-th component of the receiver to the 
    // 'icf'-th component of 'f'. 
    // A dof with coordinates (x,y,z) is initialized to f(x,y,z,t).

    main->SetTo(f,ic,icf,t) ;

    if (mortar)
	mortar->SetTo(f,ic,icf,t) ;
}


void MortaredField :: SetValues (real val)
{
    // Sets every value of the receiver to 'val'.

    main->SetValues(val) ;

    if (mortar)
	mortar->SetValues(val) ;
}


void MortaredField :: SetValues (real val, int icy)
{
    // Sets every value of the 'icy'-th component of the receiver to 'val'.

    main->SetValues(val,icy) ;

    if (mortar)
	mortar->SetValues(val,icy) ;
}

void MortaredField :: ShiftPolynomialDegree (int shift) 
{
    // Increments by 'shift' the polynomial degree of every elementary field of
    // the receiver and rebuilds the interface field accordingly.
    // The current values are lost and reinitialized to 0.
 
    main->ShiftPolynomialDegree(shift); 

    delete mortar ;
    CreateMortar() ;
}

//_____________________________ ElementaryField _______________________________


Vector<ElementaryField*>* ElementaryField :: workFields = new
Vector<ElementaryField*>(0) ;
// Initializes the static variable 'workFields'.


ElementaryField :: ElementaryField (int ncomp, boolean alloc,
                                    Vector<ParentElement*>* parents)
{
    // Constructor. Initializes the receiver to an elementary field with 'ncomp'
    // components and the parent elements specified in 'parents'.
    // Argument 'parents' is copied.
    //
    // If 'alloc' is:
    // - 'true', the values are allocated, and initialized to 0;
    // - 'false', the values are not allocated.
    //
    // On a parallel implementation, for economy of memory, only 1 processor
    // allocates the values, whereas all other data (parent elements, number of
    // components, etc) are allocated by all processors, so 'alloc' is 'true' for
    // one processor and 'false' for the others.
    // On a scalar implementation, 'alloc' should always be 'true'.

# ifdef REQUIRE
    Require("Valid number of components", ncomp >= 1) ;
    InMethod("ElementaryField::ElementaryField(ncomp,parents)") ;
# endif

    parentElements = parents->Duplicate() ;
    nbComponents   = ncomp ;
    components     = NULL ;

    if (alloc)
	AllocateValues () ;

# ifdef CHECK
    isWork = false ;
# endif
}


ElementaryField :: ElementaryField (int ncomp, boolean alloc,
                                    ParentElement* parent1)
{
    // Constructor. Initializes the receiver to an elementary field with 'ncomp'
    // components and 1 parent element ('parent1'). Such an elementary field is
    // adequate for Line elements.
    // The values are initialized to 0.
    // On the meaning of 'alloc', see the above constructor.

# ifdef REQUIRE
    Require("Valid number of components", ncomp >= 1) ;
    InMethod("ElementaryField::ElementaryField(ncomp,parent1)") ;
# endif

    parentElements = new Vector<ParentElement*>(1) ;
    parentElements->At(1) = parent1 ; 

    nbComponents = ncomp ;
    components   = NULL ;

    if (alloc)
	AllocateValues () ;

# ifdef CHECK
    isWork = false ;
# endif
}


ElementaryField :: ElementaryField (int ncomp, boolean alloc,
                                    ParentElement* parent1,
                                    ParentElement* parent2)
{
    // Constructor. Initializes the receiver to an elementary field with 'ncomp'
    // components and the parent elements 'parent1' and 'parent2'. Such an 
    // elementary field is adequate for Face elements.
    // The values are initialized to 0.
    // On the meaning of 'alloc', see the above constructor.

# ifdef REQUIRE
    Require("Valid number of components", ncomp >= 1) ;
    InMethod("ElementaryField::ElementaryField(ncomp,parent1,parent2)") ;
# endif

    parentElements = new Vector<ParentElement*>(2) ;
    parentElements->At(1) = parent1 ; 
    parentElements->At(2) = parent2 ; 

    nbComponents = ncomp ;
    components   = NULL ;

    if (alloc)
	AllocateValues () ;

# ifdef CHECK
    isWork = false ;
# endif
}


ElementaryField :: ~ElementaryField ()
{
    // Destructor.
  
    int i ;

    delete parentElements ;
    if (components) {
	for (i=1 ; i<=nbComponents ; i++)
	    delete components->At(i) ;
	delete components ;
    }
}


void ElementaryField :: Add (real alpha, int icy)
{
    // y(icy) = y(icy) + alpha .
    // Adds 'alpha' to every value of the 'icy'-th component of the receiver.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
# endif

    RealVector *val ;
    int        i, n ;

    val = components->At(icy) ;
    n   = val->GetSize() ;
    for (i=1 ; i<=n ; i++)
	val->At(i) += alpha ;
}


void ElementaryField :: AllocateValues ()
{
    // Allocates the values of the receiver. Previous values of the components,
    // if any, are lost. The values are initialized to 0.
    //
    // Dynamic allocation of the components is avoided if possible.
    //
    // On a parallel implementation, this operation is performed by 1 processor
    // only.

# ifdef REQUIRE
    Require("number of components not modified",
	    components == NULL || components->GetSize() == nbComponents) ;
    InMethod("ElementaryField::AllocateValues()") ;
# endif

    RealVector *val ;
    int        i, nval ;

    nval = 1 ;                                      // nb of values per component
    for (i=1 ; i<=parentElements->GetSize() ; i++)
	nval *= parentElements->At(i)->GetNbCollocationPoints() ;

    if (! components) {
	components = new Vector<RealVector*>(nbComponents) ;
	components->SetValues(NULL) ;
    }
    
    for (i=1 ; i<=nbComponents ; i++) {
	val = components->At(i) ;
	if (! val)
	    val = new RealVector(nval) ;
	else if (nval < val->GetSize())
	    val->SetSize(nval) ;
	else if (nval > val->GetSize())
	    val->Resize(nval) ;
	val->SetValues(ZERO) ;
	components->At(i) = val ;
    }
}


void ElementaryField :: CopyFrom (ElementaryField* x, int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y'.

# ifdef REQUIRE
    Require("y and x have same interpolation",HasSameInterpolationAs(x)) ;
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::CopyFrom(x,icx,icy)") ;
# endif

    components->At(icy)->CopyFrom(x->components->At(icx)) ;
}


void ElementaryField :: CopyInterpolateFrom (ElementaryField* x, 
                                             int icx, int icy)
{
    // y(icy) = x(icx).
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'y'.
    // If necessary, re-interpolation from 'x' to 'y' is provided.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::CopyInterpolateFrom(x,icx,icy)") ;
# endif

    TensorMatrix *interp ;
    int          i, dim ;

    if (HasSameInterpolationAs(x)) 
	components->At(icy)->CopyFrom(x->components->At(icx)) ;

    else {      // pointwise re-interpolation
	dim    = parentElements->GetSize() ;
	interp = new TensorMatrix(dim) ;
	for (i=1 ; i<=dim ; i++)
	    interp->SetMatrix(i,x->parentElements->At(i)->GetInterpMatrixTo
			      (parentElements->At(i))) ;
	interp->MatVec(x->components->At(icx),components->At(icy)) ;
	delete interp ;
    }
}
    

void ElementaryField :: CopyInterpolateTFrom (ElementaryField* x, 
                                              int icx, int icy)
{
    // y(icy) = x(icx) .
    // Copies the 'icx'-th component of 'x' into the 'icy'-th component of the
    // receiver 'x'.
    // If necessary, transposed pointwise re-interpolation from 'y' to 'x' is 
    // provided.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::CopyInterpolateTFrom(x,icx,icy)") ;
# endif

    TensorMatrix *interp ;
    int          i, dim ;

    if (HasSameInterpolationAs(x))
	components->At(icy)->CopyFrom(x->components->At(icx)) ;

    else {      // pointwise re-interpolation
	dim    = parentElements->GetSize() ;
	interp = new TensorMatrix(dim) ;
	for (i=1 ; i<=dim ; i++)
	    interp->SetMatrix(i,parentElements->At(i)->GetInterpMatrixTo
			      (x->parentElements->At(i))) ;
	interp->MatTransVec(x->components->At(icx),components->At(icy)) ;
	delete interp ;
    }
}
    

ElementaryField* ElementaryField :: CreateWork ()
{
    // Adds a new elementary field in the list of work elementary fields. 
    //
    // For generality purposes, this new elementary field is of dimension 3D; its
    // actual dimension can be modified at any time during execution.
    // Also for generality purposes, this new elementary field is given big
    // parent elements (and so a large array of values); its parent elements can
    // be modified at any time during execution.
    //
    // This method is activated when the list of work elementary fields is 
    // exhausted (i.e., contains only NULL values).

    ElementaryField        *w ;
    ParentElement          *parent ;
    Vector<ParentElement*> *vp ;
    boolean                alloc ;
    const int              BIG_SIZE = 25 ;

    parent    = ParentElementGLL::GetParentElementOfDegree(BIG_SIZE) ;
    vp        = new Vector<ParentElement*>(3) ;
    vp->At(1) = parent ;
    vp->At(2) = parent ;
    vp->At(3) = parent ;
    alloc     = true ;

    w = new ElementaryField(1,alloc,vp) ;
    workFields->Put(w) ;
    workFields->At(workFields->GetSize()) = NULL ;     // the slot of w
    delete vp ;

# ifdef CHECK
    w->isWork = true ;
# endif

    return w ;
}


void ElementaryField :: DeallocateValues ()
{
    // Deallocates the values of the receiver.
    // Used to save memory storage.
    // On a parallel implementation, this operation should be  performed by 1
    // processor only.

    int i ;

    if (components)
	for (i=1 ; i<=nbComponents ; i++) {
	    delete components->At(i) ;
	    components->At(i) = NULL ;
	}
}


void ElementaryField :: Divide (ElementaryField* x, int icx, int icy)
{
    // y = y / x .
    // Divides every value of the 'icy'-th component of the receiver 'y' by its
    // vis-a-vis in the 'icx'-th component of 'x'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::Divide(x,icx,icy)") ;
# endif

    components->At(icy)->Divide(x->components->At(icx)) ;
}


void ElementaryField :: DivideByWeights (int icy)
{
    // Divides every value of the 'icy'-th component of the receiver by the 
    // weights of the quadrature rule implemented by its parent elements.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::DivideByWeights(icy)") ;
# endif

    RealVector *w1, *w2, *w3, *val ;
    real       wj, wk, wjk ;
    int        i, j, k, dim, np1, np2, np3, idx ;
  
    dim = parentElements->GetSize() ;
    val = components->At(icy) ;

    if (dim == 1) {
	w1  = parentElements->At(1)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	for (i=1 ; i<=np1 ; i++)
	    val->At(i) /= w1->At(i) ;
    }

    else if (dim == 2) {
	w1  = parentElements->At(1)->GetWeights() ;
	w2  = parentElements->At(2)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	np2 = parentElements->At(2)->GetNbCollocationPoints() ;
	idx = 1 ;
	for (j=1 ; j<=np2 ; j++) {
	    wj = w2->At(j) ;
	    for (i=1 ; i<=np1 ; i++) {
		val->At(idx) /= w1->At(i) * wj ;
		idx++ ;
	    }
	}
    }
 
    else if (dim == 3) {
	w1  = parentElements->At(1)->GetWeights() ;
	w2  = parentElements->At(2)->GetWeights() ;
	w3  = parentElements->At(3)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	np2 = parentElements->At(2)->GetNbCollocationPoints() ;
	np3 = parentElements->At(3)->GetNbCollocationPoints() ;
	idx = 1 ;
	for (k=1 ; k<=np3 ; k++) {
	    wk = w3->At(k) ;
	    for (j=1 ; j<=np2 ; j++) {
		wjk = w2->At(j) * wk ;
		for (i=1 ; i<=np1 ; i++) {
		    val->At(idx) /= w1->At(i) * wjk ;
		    idx++ ;
		}
	    }
	}
    }

    else
	Error("ElementaryField::DivideByWeights(icy)","dim should be in {1,2,3}") ;
}


real ElementaryField :: Dot (ElementaryField* x)
{
    // dot = y(transp) . x
    // Returns the scalar product of the receiver 'y' and 'x'.
    // All components are affected.

# ifdef REQUIRE
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", nbComponents == x->nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::Dot(x)") ;
# endif

    real answer ;
    int  i ;

    answer = ZERO ;
    for (i=1 ; i<=nbComponents ; i++)
	answer += components->At(i)->Dot(x->components->At(i)) ;

    return answer ;
}


ElementaryField* ElementaryField :: Duplicate ()
{
    // Returns a new elementary field, the exact copy of the receiver.

    ElementaryField *answer ;
    int             ic ;
    boolean         alloc ;

    alloc  = (components != NULL) ;
    answer = new ElementaryField(nbComponents,alloc,parentElements) ;
    for (ic=1 ; ic<=nbComponents ; ic++)
	answer->CopyFrom(this,ic,ic) ;

    return answer ;
}


ElementaryField* ElementaryField :: DuplicateEmpty (int ncomp)
{
    // Returns a new elementary field, the exact copy of the receiver, but with
    // 'ncomp' components and all values set to 0.
    // On a parallel implementation, the returned elementary field allocates its
    // values if the receiver also does.

    ElementaryField *answer ;
    boolean         alloc ;

    if (ncomp == ALL)
	ncomp = nbComponents ;

    alloc  = (components != NULL) ;

    answer = new ElementaryField(ncomp,alloc,parentElements) ;

    return answer ;
}

real ElementaryField :: GetInfiniteNorm ()
{
    // Returns the discrete L-infinite norm (max norm) of the receiver.
    // Involves all components.

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::GetInfiniteNorm()") ;
# endif

    real answer ;
    int  i ;

    answer = ZERO ;
    for (i=1 ; i<=nbComponents ; i++)
	answer = max(answer,components->At(i)->GetInfiniteNorm()) ;

    return answer ;
}


real ElementaryField :: GetInteriorInfiniteNorm ()
{
    // Returns the L-infinite norm of the receiver. Only the internal dofs are
    // taken into account.
    // All components are involved.
    // Correct only for interpolation rules going until borders (eg, Lobatto)!

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::GetInteriorInfiniteNorm()") ;
# endif

    int dim ;
  
    dim = parentElements->GetSize() ;

    if (dim == 0)
	return GetInfiniteNorm() ;
    else if (dim == 1)
	return GetInteriorInfiniteNorm1D() ;
    else if (dim == 2)
	return GetInteriorInfiniteNorm2D() ;
    else if (dim == 3)
	return GetInteriorInfiniteNorm3D() ;
    else {
	Error("ElementaryField::GetInteriorInfiniteNorm","unvalid dimension") ;
	return 0. ;
    }
}


real ElementaryField :: GetInteriorInfiniteNorm1D ()
{
    // = method GetInteriorInfiniteNorm(), 1D case.

    RealVector *val ;
    real       answer ;
    int        i, icomp, np ;

    answer = ZERO ;
    np     = parentElements->At(1)->GetNbCollocationPoints() ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	val = components->At(icomp) ;
	for (i=2 ; i<np ; i++)
	    answer = max(answer,fabs(val->At(i))) ;
    }

    return answer ;
}


real ElementaryField :: GetInteriorInfiniteNorm2D ()
{
    // = method GetInteriorInfiniteNorm(), 2D case.

    RealVector *val ;
    real       answer ;
    int        i, j, np1, np2, nip1, icomp, idx ;

    np1    = parentElements->At(1)->GetNbCollocationPoints() ;
    np2    = parentElements->At(2)->GetNbCollocationPoints() ;
    nip1   = np1 - 2 ;
    answer = ZERO ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	val = components->At(icomp) ;
	for (j=2 ; j<np2 ; j++) {
	    idx = (j-1)*np1 + 1 ;                     // valid only if Lobatto!
	    for (i=1 ; i<=nip1 ; i++)
		answer = max(answer,fabs(val->At(idx+i))) ;
	}
    }

    return answer ;
}

real ElementaryField :: GetInteriorInfiniteNorm3D ()
{
    // = method GetInteriorInfiniteNorm(), 3D case.

    RealVector *val ;
    real       answer ;
    int        i, j, k, np1, np2, np3, nip1, icomp, idx ;

    np1    = parentElements->At(1)->GetNbCollocationPoints() ;
    np2    = parentElements->At(2)->GetNbCollocationPoints() ;
    np3    = parentElements->At(3)->GetNbCollocationPoints() ;
    nip1   = np1 - 2 ;
    answer = ZERO ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	val = components->At(icomp) ;
	for (k=2 ; k<np3 ; k++)
	    for (j=2 ; j<np2 ; j++) {
		idx = (k-1)*np1*np2 + (j-1)*np1 + 1 ;        // valid only if Lobatto!
		for (i=1 ; i<=nip1 ; i++)
		    answer = max(answer,fabs(val->At(idx+i))) ;
	    }
    }

    return answer ;
}


int ElementaryField :: GetNbDofs ()
{
    // Returns the number of collocation points per component (ie, the number of
    // values per component).

    int dir, answer ;

    answer = 1 ;
    for (dir=1 ; dir<=parentElements->GetSize() ; dir++)
	answer *= parentElements->At(dir)->GetNbCollocationPoints() ;

    return answer ;
}


int ElementaryField :: GetNbInteriorDofs ()
{
    // Returns the number of interior collocation points per component.

    int dir, answer ;

    answer = 1 ;
    for (dir=1 ; dir<=parentElements->GetSize() ; dir++)
	answer *= parentElements->At(dir)->GetNbInteriorCollocationPoints() ;

    return answer ;
}


ElementaryField* ElementaryField :: GetWork ()
{
    // Returns an elementary field with one component and with the same
    // interpolation as the receiver.

    ElementaryField *w ;
    int             i, n ;

    // 1. Find a work field

    w = NULL ;
    n = workFields->GetSize() ;
    for (i=1 ; i<=n ; i++) {
	w = workFields->At(i) ;
	if (w) {                         // found an available work field
	    workFields->At(i) = NULL ;
	    break ;
	}
    }

    if (w == NULL)                     // could not find any available work field
	w = CreateWork() ;

    // 2. Assign to w the same interpolation as that of the receiver

    n = parentElements->GetSize() ;
    w->GetParentElements()->SetSize(n) ;
    for (i=1 ; i<=n ; i++)
	w->GetParentElements()->At(i) = parentElements->At(i) ;
    w->GetComponent(1)->SetSize(components->At(1)->GetSize()) ;

    return w ;
}


ElementaryField* ElementaryField :: GetWork (ParentElement* p1, 
                                             ParentElement* p2)
{
    // Returns a 2D elementary field with one component and with parent elements
    // 'p1' and 'p2'.

    ElementaryField *w ;
    int             i, n ;

    // 1. Find a work field

    w = NULL ;
    n = workFields->GetSize() ;
    for (i=1 ; i<=n ; i++) {
	w = workFields->At(i) ;
	if (w) {                         // found an available work field
	    workFields->At(i) = NULL ;
	    break ;
	}
    }

    if (w == NULL)                     // could not find any available work field
	w = CreateWork() ;

    // 2. Assign to w the same parent elements p1 and p2

    w->GetParentElements()->SetSize(2) ;
    w->GetParentElements()->At(1) = p1 ;
    w->GetParentElements()->At(2) = p2 ;
    n = p1->GetNbCollocationPoints() * p2->GetNbCollocationPoints() ;
    w->GetComponent(1)->SetSize(n) ;

    return w ;
}


boolean ElementaryField :: HasSameInterpolationAs (ElementaryField* x)
{
    // Returns 'true' if the receiver 'y' and 'x' have the same parent elements.

    int i ;

    if (parentElements->GetSize() != x->parentElements->GetSize())
	return false ;

    else
	for (i=1 ; i<=parentElements->GetSize() ; i++)
	    if (parentElements->At(i) != x->parentElements->At(i))
		return false ;

    return true ;
}
 

real ElementaryField :: InteriorDot (ElementaryField* x)
{
    // dot = y(transp) . x
    // Returns the scalar product of the receiver 'y' and 'x'.
    // All components are affected, but only the points located in the interior
    // of the element are taken into account.
    // Correct only for interpolation rules going until borders (eg, Lobatto)!

# ifdef REQUIRE
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y and x have same nb components", nbComponents == x->nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::InteriorDot(x)") ;
# endif

    int dim ;
  
    dim = parentElements->GetSize() ;
    if (dim == 0)
	return Dot(x) ;
    else if (dim == 1)
	return InteriorDot1D(x) ;
    else if (dim == 2)
	return InteriorDot2D(x) ;
    else if (dim == 3)
	return InteriorDot3D(x) ;
    else {
	Error("ElementaryField::InteriorDot","unvalid dimension") ;
	return 0. ;
    }
}


real ElementaryField :: InteriorDot1D (ElementaryField* x)
{
    // = method InteriorDot(), 1D case.

    RealVector *valY, *valX ;
    real       answer ;
    int        icomp, np, nip ;

    answer = ZERO ;
    np     = parentElements->At(1)->GetNbCollocationPoints() ;
    nip    = np - 2 ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	valY    = components->At(icomp) ;
	valX    = x->components->At(icomp) ;
#   ifdef CHECK
	Require("correct size of valY", valY->GetSize() == np) ;
	Require("correct size of valX", valX->GetSize() == np) ;
	InMethod("ElementaryField::InteriorDot1D(x)") ;
#   endif
	answer += DDOT(&nip, &valY->At(2), &ione, &valX->At(2), &ione) ;
    }

    return answer ;
}


real ElementaryField :: InteriorDot2D (ElementaryField* x)
{
    // = method InteriorDot(), 2D case.

    RealVector *valY, *valX ;
    real       answer ;
    int        j, np1, np2, nip1, icomp, idx ;

    np1    = parentElements->At(1)->GetNbCollocationPoints() ;
    np2    = parentElements->At(2)->GetNbCollocationPoints() ;
    nip1   = np1 - 2 ;
    answer = ZERO ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	valY = components->At(icomp) ;
	valX = x->components->At(icomp) ;
#   ifdef CHECK
	Require("correct size of valY", valY->GetSize() == np1*np2) ;
	Require("correct size of valX", valX->GetSize() == np1*np2) ;
	InMethod("ElementaryField::InteriorDot2D(x)") ;
#   endif
	for (j=2 ; j<np2 ; j++) {
	    idx     = (j-1)*np1 + 2 ;
	    answer += DDOT(&nip1, &valY->At(idx), &ione, &valX->At(idx), &ione) ;
	}
    }

    return answer ;
}


real ElementaryField :: InteriorDot3D (ElementaryField* x)
{
    // = method InteriorDot(), 3D case.

    RealVector *valY, *valX ;
    real       answer ;
    int        j, k, np1, np2, np3, nip1, icomp, idx ;

    np1    = parentElements->At(1)->GetNbCollocationPoints() ;
    np2    = parentElements->At(2)->GetNbCollocationPoints() ;
    np3    = parentElements->At(3)->GetNbCollocationPoints() ;
    nip1   = np1 - 2 ;
    answer = ZERO ;
    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	valY = components->At(icomp) ;
	valX = x->components->At(icomp) ;
#   ifdef CHECK
	Require("correct size of valY", valY->GetSize() == np1*np2*np3) ;
	Require("correct size of valX", valX->GetSize() == np1*np2*np3) ;
	InMethod("ElementaryField::InteriorDot2D(x)") ;
#   endif
	for (k=2 ; k<np3 ; k++)
	    for (j=2 ; j<np2 ; j++) {
		idx     = (k-1)*np1*np2 + (j-1)*np1 + 2 ;
		answer += DDOT(&nip1, &valY->At(idx), &ione, &valX->At(idx), &ione) ;
	    }
    }

    return answer ;
}


void ElementaryField :: Inverse ()
{
    // y = 1 / y .
    // Substitutes every value of the receiver by its inverse.
    // All components are affected.

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::Inverse()") ;
# endif

    int i ;

    for (i=1 ; i<=nbComponents ; i++)
	components->At(i)->Inverse() ;
}


void ElementaryField :: Multiply (real alpha, int icy)
{
    // y = alpha*y .
    // Multiplies the 'icy'-th component of the receiver by 'alpha'.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::Multiply(alpha,icy)") ;
# endif

    components->At(icy)->Multiply(alpha) ;
}


void ElementaryField :: MultiplyAndAdd (real alpha, ElementaryField* x,
                                        int icx, int icy, real beta)
{
    // y_icy = alpha * y_icy + beta * x_icx .
    // Multiplies the 'icy'-th component of the receiver 'y' by 'alpha', then
    // adds 'beta' times the 'icx'-th component of 'x'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::MultiplyAndAdd(alpha,x,icx,icy,beta)") ;
# endif

    components->At(icy)->MultiplyAndAdd(alpha,x->components->At(icx),beta) ;
}


void ElementaryField :: MultiplyAndSwitchSigns (ElementaryField* x, 
                                                int icx, int icy)
{
    // y = - y * x .
    // Multiplies every value of the 'icy'-th component of the receiver 'y' by 
    // its vis-a-vis in the 'icx'-th component of 'x', then switches the sign of
    // 'y'. Returns 'y'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx >= 1 && icx <= x->nbComponents) ;
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("y and x have same interpolation", HasSameInterpolationAs(x)) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::MultiplyAndSwitchSigns(x,icx,icy)") ;
# endif

    components->At(icy)->MultiplyAndSwitchSigns(x->components->At(icx)) ;
}


void ElementaryField :: MultiplyByWeights (int icy)
{
    // Multiplies every value of the 'icy'-th component of the receiver by the 
    // weights of the quadrature rule implemented by its parent elements.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::MultiplyByWeights(icy)") ;
# endif

    RealVector *w1, *w2, *w3, *val ;
    real       wj, wk, wjk ;
    int        i, j, k, dim, np1, np2, np3, idx ;
  
    dim = parentElements->GetSize() ;
    val = components->At(icy) ;

    if (dim == 1) {
	w1  = parentElements->At(1)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	for (i=1 ; i<=np1 ; i++)
	    val->At(i) *= w1->At(i) ;
    }

    else if (dim == 2) {
	w1  = parentElements->At(1)->GetWeights() ;
	w2  = parentElements->At(2)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	np2 = parentElements->At(2)->GetNbCollocationPoints() ;
	idx = 1 ;
	for (j=1 ; j<=np2 ; j++) {
	    wj = w2->At(j) ;
	    for (i=1 ; i<=np1 ; i++) {
		val->At(idx) *= w1->At(i) * wj ;
		idx++ ;
	    }
	}
    }

    else if (dim == 3) {
	w1  = parentElements->At(1)->GetWeights() ;
	w2  = parentElements->At(2)->GetWeights() ;
	w3  = parentElements->At(3)->GetWeights() ;
	np1 = parentElements->At(1)->GetNbCollocationPoints() ;
	np2 = parentElements->At(2)->GetNbCollocationPoints() ;
	np3 = parentElements->At(3)->GetNbCollocationPoints() ;
	idx = 1 ;
	for (k=1 ; k<=np3 ; k++) {
	    wk = w3->At(k) ;
	    for (j=1 ; j<=np2 ; j++) {
		wjk = w2->At(j) * wk ;
		for (i=1 ; i<=np1 ; i++) {
		    val->At(idx) *= w1->At(i) * wjk ;
		    idx++ ;
		}
	    }
	}
    }

    else
	Error("ElementaryField::MultiplyByWeights(icy)","dim shoud be in {1,2,3}");
}


void ElementaryField :: OrientLike (RealVector* x)
{
    // Switches the sign of the values of the receiver 'y' wherever the 
    // condition
    //   y.x >= 0 
    // is not met.

# ifdef REQUIRE
    Require("valid dimensions", nbComponents == x->GetSize()) ;
    Require("'x' is not zero", ! x->HasOnly(ZERO)) ;
    InMethod("ElementaryField::OrientLike(x)") ;
# endif

    RealVector *y ;
    int        i, j, nbDofs ;

    y      = x->Duplicate() ;
    nbDofs = components->At(1)->GetSize() ;

    for (i=1 ; i<=nbDofs ; i++) {
	for (j=1 ; j<=nbComponents ; j++)
	    y->At(j) = components->At(j)->At(i) ;
	if (y->Dot(x) < 0)
	    // switch signs of all components of dof i
	    for (j=1 ; j<=nbComponents ; j++)
		components->At(j)->At(i) = -components->At(j)->At(i) ;
    }

    delete y ;
}


void ElementaryField :: Print ()
{
    // Prints the receiver on standard output.

    Print(true) ;
}


void ElementaryField :: Print (boolean showValues)
{
    // Prints on standard output the values of the 'icomp'-th component.

    string s ;
    int  i, dim ;

    if (components == NULL)             // not the right processor
	return ;

    if (nbComponents <= 1)
	printf("Elementary field with %d component\n",nbComponents) ;
    else
	printf("Elementary field with %d components\n",nbComponents) ;

    dim = parentElements->GetSize() ;
    if (dim) {
	printf("  parent elements are: ") ;
	for (i=1 ; i<=dim ; i++)
	  printf(" %s  ",(parentElements->At(i)->GetName(s)).c_str()) ;
	printf("\n") ;
    }

    if (showValues)
	for (i=1 ; i<=nbComponents ; i++)
	    PrintComponent(i) ;
}


void ElementaryField :: PrintComponent (int icy)
{
    // Prints on standard output the values of the 'icy'-th component.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 0 && icy <= nbComponents) ;
    InMethod("ElementaryField::PrintComponent(icy)") ;
# endif

    RealVector *comp ;
    int        i, j, k, idx, dim, n1, n2, n3 ;
    const int  MAX_PRINT = 23 ;

    if (components == NULL)                  // not the right processor
	return ;

    printf("  component %d:\n",icy) ;

    if (components->At(icy) == NULL) {       // a deleted component (saves space)
	printf("    not stored (values are 0)\n\n") ;
	return ;
    }

    comp = components->At(icy) ;
    idx  = 1  ;
    dim  = parentElements->GetSize() ;

    switch(dim) {

	case 0 :
	    printf ("%13.5E",comp->At(1)) ;
	    break ;

	case 1 :
	    n1 = parentElements->At(1)->GetNbCollocationPoints() ;
	    for (i=1; i<=n1 ; i++) {
		printf("%12.4E",comp->At(idx++)) ;
		if (i%7==0 && i!=n1)
		    printf("\n") ;
	    }
	    break ;

	case 2 :
	    n1 = parentElements->At(1)->GetNbCollocationPoints() ;
	    n2 = parentElements->At(2)->GetNbCollocationPoints() ;
	    for (j=1 ; j<=n2 ; j++) {
		//printf("    row j = %d:\n",j) ;
		for (i=1 ; i<=n1 ; i++) {
		    printf("%13.5E",comp->At(idx++)) ;
		    if (i%7==0 && i!=n1)
			printf("\n") ;
		}
		printf("\n") ;
	    }
	    break ;

	case 3 :
	    n1 = parentElements->At(1)->GetNbCollocationPoints() ;
	    n2 = parentElements->At(2)->GetNbCollocationPoints() ;
	    n3 = parentElements->At(3)->GetNbCollocationPoints() ;
	    if (n2 <= MAX_PRINT || n3 <= MAX_PRINT)
		for (k=1 ; k<=n3 ; k++) {
		    printf("    layer k = %d:\n",k) ;
		    for (j=1 ; j<=n2 ; j++) {
			printf("      row j = %d:\n",j) ;
			for (i=1 ; i<=n1 ; i++)
			    printf ("%12.4E",comp->At(idx++)) ;
			printf("\n") ;
		    }
		    printf("\n") ;
		}
	    else
		printf("     too high degree: values not printed\n\n") ;
	    break ;
    }

    printf("\n") ;
}

void ElementaryField :: Retrieve (ElementaryField* w)
{
    // Adds 'w' to the list of available work elementary fields.
    //
    // This method is a static method. So if 'w' was obtained by means of a
    // message 
    //   w = eY->GetWork() ;
    // 'w' can be retrieved by means of either of the following messages,
    // indifferently:
    //   eY->Retrieve(w) ;
    // or
    //   ElementaryField::Retrieve(w) ;

# ifdef REQUIRE
    Require ("'w' is a work field", w->isWork == true) ;
    Require ("list not full", workFields->Has(NULL)) ;
    InMethod("ElementaryField::Retrieve(w)") ;
# endif

    int i, n ;

    n = workFields->GetSize() ;
    for (i=1 ; i<=n ; i++)
	if (workFields->At(i) == NULL) {
	    workFields->At(i) = w ;
	    break ;
	}
}


void ElementaryField :: SetParentElement (ParentElement* p, int dir) 
{
    // Changes the discretization rule discretization in the 'dir'-th direction
    // to 'p'.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("has direction 'dir'", dir >= 1 && dir <= parentElements->GetSize());
    InMethod("ElementaryField::SetParentElement(p,dir)") ;
# endif

    parentElements->At(dir) = p ;

    if (components)         // do not allocate the values if the element is not
	AllocateValues() ;    //    on the current processor
}


void ElementaryField :: SetPolynomialDegree (int n)
{
    // Changes the degree of the collocation rule to 'n' in every direction.
    // The collocation rule (GL, GLL, etc) need not be the same in all
    // directions.
    // The current values are lost and reinitialized to 0.

    ParentElement *parent ;
    int           dir, ndir ;

    ndir = parentElements->GetSize() ;
    for (dir=1 ; dir<=ndir ; dir++) {
	parent = parentElements->At(dir)->GetParentElementOfDegree_(n) ;
	parentElements->At(dir) = parent ;
    }

    if (components)         // do not allocate the values if the element is not
	AllocateValues() ;    //    on the current processor
}


void ElementaryField :: SetPolynomialDegree (int n, int dir)
{
    // Changes the degree of the collocation rule to 'n' in the 'dir'-th
    // direction.
    // The collocation rule (GL, GLL, etc) need not be the same in all 
    // directions.
    // The current values are lost and reinitialized to 0.

# ifdef REQUIRE
    Require("has direction 'dir'", dir >= 1 && dir <= parentElements->GetSize());
    InMethod("ElementaryField::SetPolynomialDegree(n,dir)") ;
# endif

    ParentElement *parent ;

    parent = parentElements->At(dir)->GetParentElementOfDegree_(n) ;
    parentElements->At(dir) = parent ;

    if (components)         // do not allocate the values if the element is not
	AllocateValues() ;    //    on the current processor
}


void ElementaryField :: SetToAbsoluteValue ()
{
    // Initializes every value of the receiver to its absolute value.

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::SetToAbsoluteValue()") ;
# endif

    RealVector *comp ;
    int        i, icomp, size ;

    for (icomp=1 ; icomp<=nbComponents ; icomp++) {
	comp = components->At(icomp) ;
	size = comp->GetSize() ;
	for (i=1 ; i<=size ; i++)
	    comp->At(i) = fabs(comp->At(i)) ;
    }
}


void ElementaryField :: SetToGradient (ElementaryField* x, int icx, int icy, 
                                       int dir, ElementaryField* invJacobian)
{
    // y = dx/dX.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th coordinate (X, Y or Z)) of the 'icx'-th
    // component of 'x'.
    // For ex : in 2D : dx/dX = dx/dr dr/dX + dx/ds ds/dX.
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', interpolation from 'x' to 'y' is provided.
    // The components (dr/dX, dr/dY, ds/dX, ds/dY) of the inverse-jacobian matrix
    // are provided by 'invJacobian'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx>=1 && icx<=x->nbComponents) ;
    Require("valid argument 'icy'", icy>=1 && icy<=nbComponents) ;
    Require("valid argument 'dir'", dir>=1 && dir<=parentElements->GetSize()) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    Require("invJacobian is allocated", invJacobian->components != NULL) ;
    InMethod("FlatField::SetToGradient(x,icx,icy,dir,invJacobian)") ;
# endif

    ElementaryField *work1, *work2 ;
    int             i, dim, icInvJacobian ;

    dim   = parentElements->GetSize() ;
    work1 = GetWork() ;                             // 1 component only
    work2 = GetWork() ;                             // 1 component only

    SetValues(ZERO,icy) ;
    for (i=1 ; i<=dim ; i++) {
	icInvJacobian = dim*(i-1) + dir ;
	if (! invJacobian->HasZeroComponent(icInvJacobian)) {
	    work1->CopyInterpolateFrom(invJacobian,icInvJacobian,1) ;
	    work2->SetToGradientOnParents(x,icx,1,i) ;
	    work1->Multiply(work2,1,1) ;
	    Add(work1,1,icy) ;
	} 
    }

    Retrieve(work1) ;
    Retrieve(work2) ;
}

void ElementaryField :: SetToGradientOnParents (ElementaryField* x, 
                                                int icx, int icy, int dir)
{
    // y = dx/dr.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th local parameter (r, s or t)) of the 'icx'-th
    // component of 'x'.
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', interpolation from 'x' to 'y' is provided.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx>=1 && icx<=x->nbComponents) ;
    Require("valid argument 'icy'", icy>=1 && icy<=nbComponents) ;
    Require("valid argument 'dir'", dir>=1 && dir<=parentElements->GetSize()) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::SetToGradientOnParents(x,icx,icy,dir)") ;
# endif

    TensorMatrix *derivative ;
    int          i, dim ;

    dim        = parentElements->GetSize() ;
    derivative = new TensorMatrix(dim) ;
  
    for (i=1 ; i<=dim ; i++){
	if (x->parentElements->At(i)->GetType() != GLLI_PT){
	    if (i == dir)
		derivative->SetMatrix(i,x->parentElements->At(i)->GetInterpDerivMatrixTo
				      (parentElements->At(i))) ;
	    else
		derivative->SetMatrix(i,x->parentElements->At(i)->GetInterpMatrixTo
				      (parentElements->At(i))) ;
	}
	else{
	    if (i == dir)
		derivative->SetMatrix(i,x->parentElements->At(i)->GetDerivMatrix()) ;
	    else
		// identity matrix
		derivative->SetMatrix(i,x->parentElements->At(i)->GetInterpMatrixTo
				      (x->parentElements->At(i))) ;
	}
    }                                                 

    derivative->MatVec(x->components->At(icx),components->At(icy)) ;

    delete derivative ;
}


void ElementaryField :: SetToGradientOnParentsT (ElementaryField* x, 
                                                 int icx, int icy, int dir)
{
    // y = dx/dr.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th local parameter (r, s or t)) of the 'icx'-th
    // component of 'x'.
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', interpolation from 'x' to 'y' is provided.

# ifdef REQUIRE
    Require("valid argument 'dir'", dir>=1 && dir<=parentElements->GetSize()) ;
    Require("valid argument 'icx'", icx>=1 && icx<=x->nbComponents) ;
    Require("valid argument 'icy'", icy>=1 && icy<=nbComponents) ;
    Require("y and x have same nb components", nbComponents == x->nbComponents) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    InMethod("ElementaryField::SetToGradientOnParentsT(x,icx,icy,dir)") ;
# endif
    TensorMatrix *derivative ;
    int          i, dim ;

    dim        = parentElements->GetSize() ;
    derivative = new TensorMatrix(dim) ;

    for (i=1 ; i<=dim ; i++){
	if (x->parentElements->At(i)->GetType() != GLLI_PT){
	    if (i == dir)
		derivative->SetMatrix(i,parentElements->At(i)->GetInterpDerivMatrixTo
				      (x->parentElements->At(i))) ;
	    else
		derivative->SetMatrix(i,parentElements->At(i)->GetInterpMatrixTo
				      (x->parentElements->At(i)));
	}
	else{
	    if (i == dir)
		derivative->SetMatrix(i,parentElements->At(i)->GetDerivMatrix()) ;
	    else
		// identity matrix
		derivative->SetMatrix(i,x->parentElements->At(i)->GetInterpMatrixTo
				      (x->parentElements->At(i))) ;
	}
    } 
                                                 
    derivative->MatTransVec(x->components->At(icx),components->At(icy)) ;

    delete derivative ;

}


void ElementaryField :: SetToGradientT (ElementaryField* x, int icx, int icy,
                                        int dir, ElementaryField* invJacobian)
{
    // y = dx/dX.
    // Initializes the 'icy'-th component of the receiver 'y' to the derivative
    // (with respect to the 'dir'-th coordinate (X, Y or Z)) of the 'icx'-th
    // component of 'x'.
    // For ex : in 2D : du/dX = du/dr dr/dX + du/ds ds/dX.
    // If the icy-th parent element of 'y' differs from the icx-th parent element
    // of 'x', transposed interpolation from 'x' to 'y' is provided.
    // The components (dr/dX, dr/dY, ds/dX, ds/dY) of the inverse-jacobian matrix
    // are provided by 'invJacobian'.

# ifdef REQUIRE
    Require("valid argument 'icx'", icx>=1 && icx<=x->nbComponents) ;
    Require("valid argument 'icy'", icy>=1 && icy<=nbComponents) ;
    Require("valid argument 'dir'", dir>=1 && dir<=parentElements->GetSize()) ;
    Require("y is allocated", components != NULL) ;
    Require("x is allocated", x->components != NULL) ;
    Require("invJacobian is allocated", invJacobian->components != NULL) ;
    InMethod("FlatField::SetToGradientT(x,icx,icy,dir,invJacobian)") ;
# endif

    ElementaryField *work1, *work2 ;
    int             i, dim, icInvJacobian ;

    dim   = parentElements->GetSize() ;
    work1 = x->GetWork() ;                             // 1 component only
    work2 = GetWork() ;                                // 1 component only

    SetValues(ZERO,icy) ;
    for (i=1 ; i<=dim ; i++) {
	icInvJacobian = dim*(i-1) + dir ;
	if (! invJacobian->HasZeroComponent(icInvJacobian)) {
	    work1->CopyInterpolateFrom(invJacobian,icInvJacobian,1) ;
	    work1->Multiply(x,icx,1) ;
	    work2->SetToGradientOnParentsT(work1,1,1,i) ;
	    Add(work2,1,icy) ;
	}
    }

    Retrieve(work1) ;
    Retrieve(work2) ;
}


void ElementaryField :: SetToInvJacobian (ElementaryField* coord,
                                          ElementaryField* jacobian)
{
    // Initializes the values of the receiver to the partial derivatives of the
    // local coordinates (r,s,t) of the element with respect to the cartesian
    // coordinates (X,Y,Z).
    //
    // 'coord' and 'jacobian' provide necessary information on the element's
    // coordinates and jacobian determinant.
    //
    // See additional information in method FlatField::SetToInvJacobian().

#ifdef REQUIRE
    Require ("same interpolation as coord", HasSameInterpolationAs(coord)) ;
    Require ("same interpolation as jacobian", HasSameInterpolationAs(jacobian));
    InMethod("ElementaryField::SetToInvJacobian(coord,jacobian)") ;
# endif

    ElementaryField *grad, *work, 
	*m11, *m12, *m13, *m21, *m22, *m23, *m31, *m32, *m33 ;
    int             i, dir, dim, ic ;

    // 1. Get the partial derivatives of the coordinates

    grad = DuplicateEmpty() ;                 // better than invoking GetWork
    dim  = parentElements->GetSize() ;
    ic   = 0 ;
    for (i=1 ; i<=dim ; i++)
	for (dir=1 ; dir<=dim ; dir++) 
	    grad->SetToGradientOnParents(coord,i,++ic,dir) ;

    // 2. Invert the partial derivatives of the coordinates (later the result is
    //    multiplied by the inverse of the jacobian determinant)

    work = DuplicateEmpty(1) ;

    if (dim == 1)
	SetValues(ONE) ;

    else if (dim == 2) {                  // explicit inverse of 2-by-2 matrix
	SetValues(ZERO) ;
	Add(grad,+1.,4,1) ;                 // invJacobian(1,1) =  (2,2)
	Add(grad,-1.,2,2) ;                 // invJacobian(1,2) = -(1,2)
	Add(grad,-1.,3,3) ;                 // invJacobian(2,1) = -(2,1)
	Add(grad,+1.,1,4) ;                 // invJacobian(2,2) =  (1,1)
    }

    else if (dim == 3) {                  // explicit inverse of 3-by-3 matrix
	m11 = DuplicateEmpty(1) ;           //           by the cofactors method
	m11 ->CopyFrom(grad,5,1) ;          // m11 = (2,2)*(3,3) - (3,2)*(2,3)
	m11 ->Multiply(grad,9,1) ;
	work->CopyFrom(grad,8,1) ;
	work->Multiply(grad,6,1) ;
	m11 ->Subtract(work,1,1) ;
	m12 = DuplicateEmpty(1) ;           // m12 =-(2,1)*(3,3) + (3,1)*(2,3)
	m12 ->CopyFrom(grad,7,1) ; 
	m12 ->Multiply(grad,6,1) ;
	work->CopyFrom(grad,4,1) ;
	work->Multiply(grad,9,1) ;
	m12 ->Subtract(work,1,1) ;
	m13 = DuplicateEmpty(1) ;           // m13 = (2,1)*(3,2) - (3,1)*(2,2)
	m13 ->CopyFrom(grad,4,1) ;
	m13 ->Multiply(grad,8,1) ;
	work->CopyFrom(grad,7,1) ;
	work->Multiply(grad,5,1) ;
	m13 ->Subtract(work,1,1) ;
	m21 = DuplicateEmpty(1) ;           // m21 =-(1,2)*(3,3) + (3,2)*(1,3)
	m21 ->CopyFrom(grad,8,1) ;
	m21 ->Multiply(grad,3,1) ;
	work->CopyFrom(grad,2,1) ;
	work->Multiply(grad,9,1) ;
	m21 ->Subtract(work,1,1) ;
	m22 = DuplicateEmpty(1) ;           // m22 = (1,1)*(3,3) - (3,1)*(1,3)
	m22 ->CopyFrom(grad,1,1) ;
	m22 ->Multiply(grad,9,1) ;
	work->CopyFrom(grad,7,1) ;
	work->Multiply(grad,3,1) ;
	m22 ->Subtract(work,1,1) ;
	m23 = DuplicateEmpty(1) ;           // m23 =-(1,1)*(3,2) + (3,1)*(1,2)
	m23 ->CopyFrom(grad,7,1) ;
	m23 ->Multiply(grad,2,1) ;
	work->CopyFrom(grad,1,1) ;
	work->Multiply(grad,8,1) ;
	m23 ->Subtract(work,1,1) ;
	m31 = DuplicateEmpty(1) ;           // m31 = (1,2)*(2,3) - (2,2)*(1,3)
	m31 ->CopyFrom(grad,2,1) ;
	m31 ->Multiply(grad,6,1) ;
	work->CopyFrom(grad,5,1) ;
	work->Multiply(grad,3,1) ;
	m31 ->Subtract(work,1,1) ;
	m32 = DuplicateEmpty(1) ;           // m32 =-(1,1)*(2,3) + (2,1)*(1,3)
	m32 ->CopyFrom(grad,4,1) ;
	m32 ->Multiply(grad,3,1) ;
	work->CopyFrom(grad,1,1) ;
	work->Multiply(grad,6,1) ;
	m32 ->Subtract(work,1,1) ;
	m33 = DuplicateEmpty(1) ;           // m33 = (1,1)*(2,2) - (2,1)*(1,2)
	m33 ->CopyFrom(grad,1,1) ;
	m33 ->Multiply(grad,5,1) ;
	work->CopyFrom(grad,4,1) ;
	work->Multiply(grad,2,1) ;
	m33 ->Subtract(work,1,1) ;
	CopyFrom(m11,1,1) ;                 // invJacobian(1,1) = m11
	CopyFrom(m21,1,2) ;                 // invJacobian(1,2) = m21
	CopyFrom(m31,1,3) ;                 // invJacobian(1,3) = m31
	CopyFrom(m12,1,4) ;                 // invJacobian(2,1) = m12
	CopyFrom(m22,1,5) ;                 // invJacobian(2,2) = m22
	CopyFrom(m32,1,6) ;                 // invJacobian(2,3) = m32
	CopyFrom(m13,1,7) ;                 // invJacobian(3,1) = m13
	CopyFrom(m23,1,8) ;                 // invJacobian(3,2) = m23
	CopyFrom(m33,1,9) ;                 // invJacobian(3,3) = m33
	delete m11 ;
	delete m12 ;
	delete m13 ;
	delete m21 ;
	delete m22 ;
	delete m23 ;
	delete m31 ;
	delete m32 ;
	delete m33 ;
    }

    else
	Error("ElementaryField::SetToInvJacobian(coord,jacobian)",
	      "mesh dimension should be 1, 2 or 3") ;
 
    // 3. Divide every component by the jacobian

    work->CopyFrom(jacobian,1,1) ;                       // work = 1/jac
    work->Inverse() ;
 
    ic = 0 ;
    for (i=1 ; i<=dim ; i++)
	for (dir=1 ; dir<=dim ; dir++)
	    Multiply(work,1,++ic) ;                          // multiply by 1/jac

    delete work ;
    delete grad ;

    // 4. Something ugly that sets to 0 the components that are almost 0, in
    //    order to execution time in later calculations and also memory space

#ifdef FAST_GRADIENTS

    real       largestNorm, norm, ratio ;
    const real PREC = 1.e-12 ;

    largestNorm = ZERO ;
    for (ic=1 ; ic<=dim*dim ; ic++) 
	largestNorm=max(largestNorm,GetComponent(ic)->GetEuclidianNorm()) ;
    for (ic=1 ; ic <= dim*dim ; ic++) {
	norm  = GetComponent(ic)->GetEuclidianNorm() ;
	ratio = norm / largestNorm ;
	if (ratio < PREC)
	    // in order to save memory, component ic is not stored
	    SetZeroComponent(ic) ;
    }

# endif
}


void ElementaryField :: SetToJacobian (ElementaryField* coord)
{
    // Initializes the values of the receiver to the jacobian, i.e., to the
    // absolute value of the determinant of the matrix of the partial derivatives
    // of the cartesian coordinates (X,Y,Z) with respect to the local coordinates
    // (r,s,t) of the element.
    //
    // More specifically:
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
    //
    // 'coord' provides information on the element's coordinates.

#ifdef REQUIRE
    Require("the field has 1 component", nbComponents == 1) ;
    InMethod("ElementaryField::SetToJacobian(coord)") ;
# endif

    ElementaryField *grad, *work, *m11, *m21, *m31, *N, *jac ;
    int             i, dir, dim, ic, nbCoord ;

    dim     = parentElements->GetSize() ;
    nbCoord = coord->GetNbComponents() ;

//  printf ("dim:%d\n",dim);
//  printf ("nbCoord:%d\n",nbCoord);

    // 1. Cases 4 and 5

    if ((dim == 1 && nbCoord == 2) || (dim == 2 && nbCoord == 3)) {

	N   = DuplicateEmpty(nbCoord) ;
	jac = N->SetToNormal(coord) ;        // jac is a by-product of SetToNormal
	this->CopyFrom(jac,1,1) ;
	delete N ;
	delete jac ;
    }

    // 2. Cases 1, 2 and 3

    else {

	// 2.a. Get the partial derivatives of the coordinates

	grad = DuplicateEmpty(dim*dim) ;
	ic   = 0 ;
	for (i=1 ; i<=dim ; i++)
	    for (dir=1 ; dir<=dim ; dir++)
		grad->SetToGradientOnParents(coord,i,++ic,dir) ;

	// 2.b. Compute the determinant

	if (dim == 1 && nbCoord == 1) {          // jac = (1,1)
	    // case 1
	    CopyFrom(grad,1,1) ;
	    SetToAbsoluteValue() ;
	}

	else if (dim == 2 && nbCoord == 2) {     // jac = (1,1)*(2,2) + (2,1)*(1,2)
	    // case 2
	    work = DuplicateEmpty(1) ;
	    this->CopyFrom(grad,1,1) ;
	    this->Multiply(grad,4,1) ;
	    work->CopyFrom(grad,2,1) ;
	    work->Multiply(grad,3,1) ;
	    this->Subtract(work,1,1) ;
	    SetToAbsoluteValue() ;
	    delete work ;
	}
  
	else if (dim == 3 && nbCoord == 3) {     // jac = (1,1)*m11 - (2,1)*m21
	    // case 3                              //     + (3,1)*m31
	    work = DuplicateEmpty(1) ;

	    m11 = DuplicateEmpty(1) ;              // m11 = (2,2)*(3,3) - (3,2)*(2,3)
	    m11 ->CopyFrom(grad,5,1) ;
	    m11 ->Multiply(grad,9,1) ;
	    work->CopyFrom(grad,8,1) ;
	    work->Multiply(grad,6,1) ;
	    m11 ->Subtract(work,1,1) ;
  
	    m21 = DuplicateEmpty(1) ;              // m21 =-(1,2)*(3,3) + (3,2)*(1,3)
	    m21 ->CopyFrom(grad,8,1) ;
	    m21 ->Multiply(grad,3,1) ;
	    work->CopyFrom(grad,2,1) ;
	    work->Multiply(grad,9,1) ;
	    m21 ->Subtract(work,1,1) ;
  
	    m31 = DuplicateEmpty(1) ;              // m31 = (1,2)*(2,3) - (2,2)*(1,3)
	    m31 ->CopyFrom(grad,2,1) ;
	    m31 ->Multiply(grad,6,1) ;
	    work->CopyFrom(grad,5,1) ;
	    work->Multiply(grad,3,1) ;
	    m31 ->Subtract(work,1,1) ;
  
	    this->CopyFrom(grad,1,1) ;             // jac
	    this->Multiply(m11,1,1) ;
	    work->CopyFrom(grad,4,1) ;
	    work->Multiply(m21) ;
	    this->Add(work,1,1) ;
	    work->CopyFrom(grad,7,1) ;
	    work->Multiply(m31) ;
	    this->Add(work,1,1) ;
	    SetToAbsoluteValue() ;
  
	    delete m11 ;
	    delete m21 ;
	    delete m31 ;
	    delete work ;
	}

	else
	    Error("ElementaryField::SetToJacobian(coord)","unexpected case") ;
   
	delete grad ;
    }
}

ElementaryField* ElementaryField :: SetToNormal (ElementaryField* coord)
{
    // Initializes the values of the receiver to the normal vector "n", which is
    // a function of the partial derivatives of the cartesian coordinates (X,Y,Z)
    // with respect to the local coordinates (r,s,t) of the element.
    // At every dof, n has euclidian norm equal to 1.
    //
    // n is given by the expression N = J n, where:
    // - J is the jacobian, equal to the euclidian norm of N,
    // - N is the non-normalized normal vector given by:
    //    . case a: 1D element in a 2D space
    //       N = (-dY/dr, dX/dr) = orthogonal to (dX/dr, dY/dr)
    //    . case b: 2D element in a 3D space
    //       N = cross product of (dX/dr, dY/dr, dZ/dr) and (dX/ds, dY/ds, dZ/ds)
    //
    // Returns a new elementary field initialized to J, as a by-product of the
    // computations.
    //
    // Remark. This method does set any particular orientation to "n", i.e, the
    //         sign of the returned values is arbitrary. Therefore the calling
    //         method must include additional provisions in order to orient "n"
    //         adequately.

    ElementaryField *grad, *work, *jac ;
    int             i, dir, dim, ic, nbCoord ;

    // 1. Get the partial derivatives of the coordinates

    dim     = parentElements->GetSize() ;
    nbCoord = coord->GetNbComponents() ;

    grad    = DuplicateEmpty(nbCoord*dim) ;
    ic      = 0 ;
    for (i=1 ; i<=nbCoord ; i++)
	for (dir=1 ; dir<=dim ; dir++)
	    grad->SetToGradientOnParents(coord,i,++ic,dir) ;

    // 2. Get the non-normalized normal vector N

    jac  = DuplicateEmpty(1) ;
    work = DuplicateEmpty(1) ;

    if (dim == 1 && nbCoord == 2) {          // N = ( -(1,2), (1,1) )
	// case a 

	this->CopyFrom(grad,2,1) ;             // N(1) = -(1,2)
	this->Multiply(ONE_MINUS,1) ;

	this->CopyFrom(grad,1,2) ;             // N(2) =  (1,1)

	jac->CopyFrom(this,1,1) ;              // N(1) * N(1)
	jac->Multiply(jac) ;
	work->CopyFrom(this,2,1) ;             // N(2) * N(2)
	work->Multiply(work) ;
	jac->Add(work,1,1) ;
	jac->Sqrt() ;                          // jac  = norm(N)
    }

    else if (dim == 2 && nbCoord == 3) {     // N = ((1,1), (1,2), (1,3)) cross
	// case b                              //     ((2,1), (2,2), (2,3))

	this->CopyFrom(grad,3,1) ;             // N(1) = (1,2)*(2,3) - (1,3)*(2,2)
	this->Multiply(grad,6,1) ;
	work->CopyFrom(grad,4,1) ;
	work->Multiply(grad,5,1) ;
	this->Subtract(work,1,1) ;

	this->CopyFrom(grad,2,2) ;             // N(2) = (1,3)*(2,1) - (1,1)*(2,3)
	this->Multiply(grad,5,2) ;
	work->CopyFrom(grad,1,1) ;
	work->Multiply(grad,6,1) ;
	this->Subtract(work,1,2) ;

	this->CopyFrom(grad,1,3) ;             // N(3) = (1,1)*(2,2) - (1,2)*(2,1)
	this->Multiply(grad,4,3) ;
	work->CopyFrom(grad,2,1) ;
	work->Multiply(grad,3,1) ;
	this->Subtract(work,1,3) ;

	jac->CopyFrom(this,1,1) ;              // N(1) * N(1)
	jac->Multiply(jac) ;
	work->CopyFrom(this,2,1) ;             // N(2) * N(2)
	work->Multiply(work,1,1) ;
	jac->Add(work,1,1) ;
	work->CopyFrom(this,3,1) ;             // N(3) * N(3)
	work->Multiply(work,1,1) ;
	jac->Add(work,1,1) ;
	jac->Sqrt() ;                          // jac  = norm(N)
    }

    else
	Error("ElementaryField::SetToNormal(coord)","unexpected case") ;

    // 3. Normalization: n = N / jac

    for (i=1 ; i<=nbCoord ; i++)
	this->Divide(jac,1,i) ;                // n(i) = N(i) / jac

    delete grad ;
    delete work ;

    return jac ;
}

void ElementaryField :: SetValues (real val)
{
    // Sets every value of every component of the receiver to 'val'.

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::SetValues(val)") ;
# endif

    int i ;

    if (val == ZERO)
	for (i=1 ; i<=nbComponents ; i++)
	    components->At(i)->SetValuesZero() ;

    else
	for (i=1 ; i<=nbComponents ; i++)
	    components->At(i)->SetValues(val) ;
}


void ElementaryField :: SetValues (real val, int icy)
{
    // Sets every value of the 'icy'-th component of the receiver to 'val'.

# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::SetValues(val,icy)") ;
# endif

    if (val == ZERO)
	components->At(icy)->SetValuesZero() ;
    else
	components->At(icy)->SetValues(val) ;
}


void ElementaryField :: SetZeroComponent (int icy)
{
    // Sets to zero the values of the 'icy'-th component of the receiver.
    // Furthermore, in order to save memory space, the component is physically
    // deleted. 
    // This is useful to avoid the storage of the nul components of the inverse
    // jacobian matrix of a straight-geometry element. For example, for a
    // parallelipipedic 3D element, the storage of 6 of the 9 components can be
    // saved.
    // The values of the component should never be used in any further 
    // computations, otherwise a crash may occur!
  
# ifdef REQUIRE
    Require("valid argument 'icy'", icy >= 1 && icy <= nbComponents) ;
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::SetZeroComponent(icy)") ;
# endif

    delete components->At(icy) ;
    components->At(icy) = NULL ;
}


void ElementaryField :: ShiftPolynomialDegree (int shift) 
{
    // Increments by 'shift' the polynomial degree of the parent element in every
    // direction.
    // The collocation rule and the polynomial degree need not be identical in
    // all directions.
    // The current values are lost and reinitialized to 0.
 
    ParentElement *newParent ;
    int           dir, ndir, newDegree ;

    ndir = parentElements->GetSize() ;
    for (dir=1 ; dir<=ndir ; dir++) {
	newDegree = parentElements->At(dir)->GetPolynomialDegree() + shift ;
	newParent = parentElements->At(dir)->GetParentElementOfDegree_(newDegree) ;
	parentElements->At(dir) = newParent ;
    }

    if (components)         // do not allocate the values if the element is not
	AllocateValues() ;    //    on the current processor
}

void ElementaryField :: Square ()
{ 
    // y = square(y).
    // Replaces every value of the receiver by its square.

    int i ;
  
    for (i=1 ; i<=nbComponents ; i++)
	components->At(i)->Square() ;
}

void ElementaryField :: Sqrt ()
{
    // y = sqrt(y).
    // Replaces every value of the receiver by its square root. All values of the
    // receiver should be positive.

    int i ;
  
    for (i=1 ; i<=nbComponents ; i++)
	components->At(i)->Sqrt() ;
}

void ElementaryField :: Power(real x)
{
    // MAH 04.2006

    int i ;
  
    for (i=1 ; i<=nbComponents ; i++)
	components->At(i)->Power(x) ;
}


real ElementaryField :: SumValues ()
{
    // Returns the sum of all values of all components of the receiver.

# ifdef REQUIRE
    Require("is allocated", components != NULL) ;
    InMethod("ElementaryField::GetInfiniteNorm()") ;
# endif

    real answer ;
    int  i ;

    answer = ZERO ;
    for (i=1 ; i<=nbComponents ; i++)
	answer += components->At(i)->SumValues() ;

    return answer ;
}


void ElementaryField :: ScalarProductWithTensors (ElementaryField* x1, ElementaryField* x2)
{   // MAH 12.2005
    // Calcule le produit scalaire x_ij.y_ij et le met dans le receveur.


#ifdef REQUIRE
    Require("receveur nbcomponents=1.",nbComponents == 1) ;
    Require("x1 et x2 -> Nb de composantes egale.",x1->GetNbComponents() == x2->GetNbComponents()) ;
    InMethod("ElementaryField::ScalarProductWithTensors(x1,x2)") ;
# endif

    int i, j, nb, dim ;
    int n1=1, n2=1, n3=1 ;

    dim = parentElements->GetSize() ;

    n1 = (dim >= 1) ? parentElements->At(1)->GetNbCollocationPoints() : 1 ;
    n2 = (dim >= 2) ? parentElements->At(2)->GetNbCollocationPoints() : 1 ;
    n3 = (dim == 3) ? parentElements->At(3)->GetNbCollocationPoints() : 1 ;
  
    nb = n1*n2*n3; 

//for (j=1 ; j<=nb; j++) 
//      components->At(1)->At(j) = 0. ;

    for (j=1 ; j<=nb; j++) 
	for (i=1 ; i<=x1->GetNbComponents(); i++){
	    components->At(1)->At(j) += x1->GetComponent(i)->At(j) *
		x2->GetComponent(i)->At(j) ;
	}   
}


void ElementaryField :: SetToTrace (ElementaryField* x)
{   // MAH 07.12.2005
    // Receiver gets the trace of x
    // x must be a tensor of order 2 with 9 components

#ifdef REQUIRE
    Require("receveur nbcomponents=1.",nbComponents == 1) ;
    InMethod("ElementaryField::SetToTrace(x)") ;
# endif

    int i, j, k,  nb, dim ;
    int n1=1, n2=1, n3=1 ;

    dim = parentElements->GetSize() ;

    n1 = (dim >= 1) ? parentElements->At(1)->GetNbCollocationPoints() : 1 ;
    n2 = (dim >= 2) ? parentElements->At(2)->GetNbCollocationPoints() : 1 ;
    n3 = (dim == 3) ? parentElements->At(3)->GetNbCollocationPoints() : 1 ;
  
    nb = n1*n2*n3; 

    for (j=1 ; j<=nb; j++) 
	for (i=1 ; i<=3; i++){
	    k=(i-1)*3+i ;
	    components->At(1)->At(j) += x->GetComponent(k)->At(j) ;
	}   
}
