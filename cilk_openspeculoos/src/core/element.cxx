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

// element.cxx

#include "core/element.hxx"
#include "core/point.hxx"
#include "core/vertex.hxx"
#include "core/edge.hxx"
#include "core/face.hxx"
#include "core/volume.hxx"
#include "core/mesh.hxx"
#include "core/distrib.hxx"
#include "core/matrix.hxx"
#include "core/funct.hxx"
#include "core/util.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <stdio.h>


Element :: Element ()
    : GarbagedObject()
{
    // Constructor.

    distributions    = NULL ;
    elementaryFields = new Vector<ElementaryField*>(0) ; 
# ifdef PARALLEL
    procNumber       = UNINITIALIZED_INTEGER ;
# else
    procNumber       = 0 ;
# endif
}


Element :: ~Element () 
{
    // Destructor.

    int i ;

    for (i=1 ; i<=distributions->GetSize() ; i++)
	unref(distributions->At(i)) ;
    delete distributions ;

    for (i=1 ; i<=elementaryFields->GetSize() ; i++)
	delete elementaryFields->At(i) ;
    delete elementaryFields ;
}


void Element :: CopyAddValuesFrom (Element* elemX, int indX, int indY,
                                   int icx, int icy, CopyAddMode cam, 
                                   InterpMode im, NonConform nc,
                                   ReceiveMode rm)
{
    // Copies or adds values from the 'indX'-th elementary field of 'elemX' to
    // the 'indY'-th elementary field of the receiver.
    // See method Field::CopyAddValuesFrom(x,icx,icy,cam,im,nc,ct).

# ifdef REQUIRE
    Require("valid argument 'indY'", HasElementaryField(indY)) ;
    Require("valid argument 'indX'", elemX->HasElementaryField(indX)) ;
    Require("consistent arguments icx and icy",
	    (icx==ALL && icy==ALL) || (icx!=ALL && icy!=ALL)) ;
    InMethod("Element::CopyAddValuesFrom(elemX,indX,indY,icx,icy,cam,im,nc)") ;
# endif
 
    ElementaryField *eY, *eX ;
    int             i, dim, nbComp ;

    eY  = GetElementaryField(indY) ;
    eX  = elemX->GetElementaryField(indX) ;
    dim = elemX->GetDimension() ;

    if (icx == ALL && icy == ALL) {
	nbComp = eY->GetNbComponents() ;
#   ifdef CHECK
	if (eX->GetNbComponents() != eY->GetNbComponents())
	    Error("Element::CopyAddValuesFrom","ncx != ncy : unexpected case") ;
#   endif
	for (i=1 ; i<=nbComp ; i++) {
	    if (dim == 0)
		CopyAddValuesFromVertex((Vertex*)elemX,indX,indY,i,i,cam,im,nc,rm) ;
	    else if (dim == 1)
		CopyAddValuesFromEdge((Edge*)elemX,indX,indY,i,i,cam,im,nc,rm) ;
	    else if (dim == 2)
		CopyAddValuesFromFace((Face*)elemX,indX,indY,i,i,cam,im,nc,rm) ;
	    else
		CopyAddValuesFromVolume((Volume*)elemX,indX,indY,i,i,cam,im,nc,rm) ;
	}
    }

    else {
	if (dim == 0)
	    CopyAddValuesFromVertex((Vertex*)elemX,indX,indY,icx,icy,cam,im,nc,rm) ;
	else if (dim == 1)
	    CopyAddValuesFromEdge((Edge*)elemX,indX,indY,icx,icy,cam,im,nc,rm) ;
	else if (dim == 2) 
	    CopyAddValuesFromFace((Face*)elemX,indX,indY,icx,icy,cam,im,nc,rm) ;
	else
	    CopyAddValuesFromVolume((Volume*)elemX,indX,indY,icx,icy,cam,im,nc,rm) ;
    }
}


void Element :: CopyAddValuesFromEdge (Edge*, int, int, int, int, 
                                       CopyAddMode, InterpMode, NonConform,
                                       ReceiveMode)
{
    // Issues an error message.

    ImplementedBySubclasses("Element::CopyAddValuesFromEdge") ;
}


void Element :: CopyAddValuesFromFace (Face*, int, int, int, int, 
                                       CopyAddMode, InterpMode, NonConform,
                                       ReceiveMode)
{
    // Issues an error message.

    ImplementedBySubclasses("Element::CopyAddValuesFromFace") ;
}


void Element :: CopyAddValuesFromVertex (Vertex*, int, int, int, int, 
                                         CopyAddMode, InterpMode, NonConform,
                                         ReceiveMode)
{
    // Issues an error message.

    ImplementedBySubclasses("Element::CopyAddValuesFromVertex") ;
}


void Element :: CopyAddValuesFromVolume (Volume*, int, int, int, int, 
                                         CopyAddMode, InterpMode, NonConform,
                                         ReceiveMode)
{
    // Issues an error message.

    ImplementedBySubclasses("Element::CopyAddValuesFromVolume") ;
}


void Element :: CreateElementaryField (int index, int ncomp, 
                                       Vector<ParentElement*>* vp)
{
    // Adds to the list of elementary fields of the receiver a new elementary
    // field with index 'index', 'ncomp' components and a copy of 'vp' as list of
    // parent elements.
    // The size of 'vp' is usually equal to the geometrical dimension of the
    // receiver; in rare cases, e.g., coordinate field of a skeleton mesh, the
    // size of 'vp' is larger.

# ifdef REQUIRE
    Require("new elementary field", ! HasElementaryField(index)) ; 
    Require("valid number of components", ncomp > 0) ;
    InMethod("Element::CreateElementaryField(index,ncomp,vp)") ;
# endif

    boolean alloc ;

    if (elementaryFields->GetSize() < index)
	elementaryFields->Resize(index) ;

    alloc = (Mpi_rank == procNumber) ;
    elementaryFields->At(index) = new ElementaryField(ncomp,alloc,vp) ;
}


void Element :: DeleteElementaryField (int index)
{
    // Deletes the elementary field with index 'index' from the list of 
    // elementary fields of the receiver.

# ifdef REQUIRE
    Require("elementary field exists", HasElementaryField(index)) ; 
    InMethod("Element::DeleteElementaryField(index)") ;
# endif

    delete elementaryFields->At(index) ;
    elementaryFields->At(index) = NULL ;
}


void Element :: DuplicateElementaryField (int index1, int index2)
{
    // Creates and stores a new elementary field with index 'index2' as a copy 
    // of the existing elementary field with index 'index1'.

# ifdef REQUIRE
    Require("elem. field index1 exists", HasElementaryField(index1)) ; 
    Require("elem. field index2 does not exist", ! HasElementaryField(index2)) ; 
    InMethod("Element::DuplicateElementaryField(index1,index2)") ;
# endif

    if (elementaryFields->GetSize() < index2)
	elementaryFields->Resize(index2) ;

    elementaryFields->At(index2) = elementaryFields->At(index1)->Duplicate() ;
}


void Element :: DuplicateEmptyElementaryField (int index1, int index2,
                                               int ncomp)
{
    // Creates and stores a new elementary field with index 'index2' as a copy 
    // of the existing elementary field with index 'index1', but with 'ncomp'
    // components and all values initialized to 0.

# ifdef REQUIRE
    Require("elem. field index1 exists", HasElementaryField(index1)) ; 
    Require("elem. field index2 does not exist", ! HasElementaryField(index2)) ; 
    InMethod("Element::DuplicateElementaryField(index1,index2,ncomp)") ;
# endif

    if (elementaryFields->GetSize() < index2)
	elementaryFields->Resize(index2) ;

    elementaryFields->At(index2) = elementaryFields->At(index1)
	->DuplicateEmpty(ncomp) ;
}


void Element:: PrintBarycenter ()
{ 
  // Prints on standard output the barycenter of every element in the receiver.
  // Less verbose than method 'Print'.

    Point   *bary ;

    if (Mpi_rank == 0) {
	printf("\nBarycenter de l'element mortier:\n");
	bary = GetBarycenter() ;
	bary->Print() ;
	delete bary ;
    }
}

Vector<ParentElement*>* Element :: GetCommonParentElements (Element* elemX,
                                                            int indX)
{
    // The receiver is assumed to be included in 'elemX'. Among the parent 
    // elements of the 'indx'-elementary field of 'elemX', returns those which
    // correspond to the same directions as the receiver.
    // Programming note: this method could be derived into classes Edge, Face
    //    and Volume. This is not done here only for keeping it compact.

# ifdef REQUIRE
    Require("'elemX' contains the receiver", this->IsIncludedIn(elemX)) ;
    Require("valid dimensions", elemX->GetDimension() == GetDimension() ||
	    elemX->GetDimension() == GetDimension()+1) ;
    InMethod("Element::GetCommonParentElements(elemX,indx,indy)") ;
# endif

    Vector<ParentElement*> *answer ;
    ElementaryField        *eX ;
    int                    dir, dir1, dir2, dimX ;

    dimX = elemX->GetDimension() ;
    eX   = elemX->GetElementaryField(indX) ;

    if (dimX == GetDimension())     // most probably, receiver = elemX
	answer = eX->GetParentElements()->Duplicate() ;

    else {
	answer = new Vector<ParentElement*>(GetDimension()) ;
	if (dimX == 1) {
	    ;                           // do nothing. VVK wrote: ParentElementGLL(0)
	}
	else if (dimX == 2) {
	    dir = ((Face*)elemX)->GetDirectionOf((Edge*)this) ;
	    answer->At(1) = eX->GetParentElement(dir) ;
	}
	else if (dimX == 3) {
	    ((Volume*)elemX)->GetDirectionsOf((Face*)this,dir1,dir2) ;
	    answer->At(1) = eX->GetParentElement(dir1) ;
	    answer->At(2) = eX->GetParentElement(dir2) ;
	}
    }

    return answer ;
}


int Element :: GetDimension ()
{
    // Strangely enough, cannot be declared as pure virtual, since it is called
    // in the destructor.

    Error("Element::GetDimension","pure virtual method") ;
    return 0 ;
}


Distribution* Element :: GetDistribution (int i)
{
    // Returns the refinement of the receiver in the 'i'-th direction. If such 
    // distribution does not exist yet, creates it as a uniform distribution of 
    // size 1.

# ifdef REQUIRE
    Require("valid index 'i'", i>=0 && i<=GetDimension()) ;
    InMethod("Element::GetDistribution(i)") ;
# endif

    if (distributions->At(i) == NULL)
	SetDistribution(i,1) ;

    return distributions->At(i) ;
}


int Element :: GetNbSubelements ()
{
    // Issues an error message.

    ImplementedBySubclasses("Element::GetNbSubelements") ;
    return 0 ;
}


void Element :: GetSetDof (int index, int iop, int idof, real &val)
{
    // Gets ('iop'=1) or sets ('iop'=0) the 'idof'-th degree of freedom to 'val',
    // for the 'index'-th field of the receiver.
    // The numbering follows the rule: starting from 1, the indices run faster in
    // direction 1, then in the subsequent direction and finally on the number of
    // components.

    RealVector *values ;
    int        size, ic, ival ;

    size   = GetElementaryField(index)->GetNbDofs() ;
    ic     = (idof-1)/size + 1 ;
    ival   = (idof-1)%size + 1 ;
    values = GetElementaryField(index)->GetComponent(ic) ;

    if (iop)
	values->At(ival) = val ;
    else
	val = values->At(ival) ;
}


void Element :: GetSetInteriorDof (int, int, int, real&)
{
    // Implemented by subclasses.

    ImplementedBySubclasses("Element::GetSetInteriorDof") ;
}


Element* Element :: GetSubelement (int)
{
    // Implemented by subclasses.

    ImplementedBySubclasses("Element::GetSubelement(i)") ;
    return NULL ;
}


real Element :: GetTypicalLength (int)
{
    // Issues an error message.

    ImplementedBySubclasses("Element::GetTypicalLength(dir)") ;
    return 0. ;
}


boolean Element :: Has (Element*)
{
    // Implemented by subclasses.

    ImplementedBySubclasses("Element::Has(elem)") ;
    return 0 ;
}


void Element :: InitializeDistributions ()
{
    // Creates attribute "distributions".
    // C++ is a shit: this method should be part of the constructor of this
    // class. However, message GetDimension is not executed with late binding:
    // the dummy version of class Element is invoked. why the hell?

    int i, n ;

    n             = GetDimension() ;
    distributions = new Vector<Distribution*>(n) ;
    for (i=1 ; i<=n ; i++)
	distributions->At(i) = NULL ;
}


boolean Element :: IsEdge ()
{
    // Returns 'false', because the receiver is not an edge.

    return false ;
}


boolean Element :: IsFace ()
{
    // Returns 'false', because the receiver is not a face.

    return false ;
}


boolean Element :: IsIdenticalTo (Element*)
{
    // Implemented by subclasses.

    ImplementedBySubclasses("Element::IsIdenticalTo(elem)") ;
    return false ;
}


boolean Element :: IsVertex ()
{
    // Returns 'false', because the receiver is not a vertex.

    return false ;
}


boolean Element :: IsVolume ()
{
    // Returns 'false', because the receiver is not a volume.

    return false ;
}


void Element :: PrintPointers ()
{
    // Implemented by subclasses.

    ImplementedBySubclasses("Element::PrintPointers()") ;
}


void Element :: SetDistribution (int i, Distribution* d)
{
    // Sets the refinement of the receiver in the i-th direction to 'd'.
    // Does not make the receiver owner of 'd'.

# ifdef REQUIRE
    Require("valid argument 'i'", i>=1 && i<=GetDimension()) ;
    InMethod("Element::SetDistribution(i,d)") ;
# endif

    ref(d) ;
    unref(distributions->At(i)) ;
    distributions->At(i) = d ;
}


void Element :: SetDistribution (int i, int n)
{
    // Sets the refinement of the receiver in the i-th direction to a uniform 
    // distribution of size 'n'.

# ifdef REQUIRE
    Require("valid argument 'i'", i>=1 && i<=GetDimension()) ;
    Require("valid argument 'n'", n>=1) ;
    InMethod("Element::SetDistribution(i,d)") ;
# endif

    unref(distributions->At(i)) ;
    distributions->At(i) = new UniformDistribution(n) ;
}


void Element :: SetProcessorNumber (int p)
{
    // Assigns the receiver to the processor of rank 'p'.
    // This can be done only before the receiver has created any elementary
    // field.
    // On a 1-processor implementation, invoking this method is either useless
    // (if 'p' is 0) or harmful (if 'p' is not 0).

# ifdef REQUIRE
    Require("has no elementary fields yet", elementaryFields->GetSize() == 0) ;
    InMethod("Element::SetProcessorNumber(p)") ;
# endif

    procNumber = p ;
}


void Element :: SetTo (int index, int indexCoord, Function* f, int ic, int icf,
                       real t)
{
    // Initializes the 'ic'-th component of the 'index'-th elementary field of
    // the receiver to the values of the 'icf'-th component of 'f'. A dof with
    // coordinates (x,y,z,t) is initialized to f(x,y,z,t).
    // Argument 'indexCoord' gives the index of the coordinate field (x,y,z);
    // this field is expected to have the same interpolation as the receiver.

    ElementaryField *coord ;
    RealVector      *val, *valX, *valY, *valZ ;
    int             i, ncoord, size ;

# ifdef REQUIRE
    Require("valid argument 'ic'", 
	    ic>=1 && ic<=GetElementaryField(index)->GetNbComponents()) ;
    Require("field and coordinate field have same interpolation",
	    GetElementaryField(indexCoord)->HasSameInterpolationAs
	    (GetElementaryField(index))) ;
    InMethod("Element::SetTo(index,indexCoord,f,ic,icf,t)") ;
# endif

    if (Mpi_rank != procNumber)
	return ;

    coord  = GetElementaryField(indexCoord) ;
    ncoord = coord->GetNbComponents() ;
    val    = GetElementaryField(index)->GetComponent(ic) ;
    valX   = coord->GetComponent(1) ;
    size   = valX->GetSize() ;

    if (ncoord == 1) {
	for (i=1 ; i<=size ; i++) 
	    val->At(i) = f->Evaluate(valX->At(i),ZERO,ZERO,t,icf) ;
    }

    else if (ncoord == 2) {
	valY = coord->GetComponent(2) ;
	for (i=1 ; i<=size ; i++)
	    val->At(i) = f->Evaluate(valX->At(i),valY->At(i),ZERO,t,icf) ;
    }

    else if (ncoord == 3) {
	valY = coord->GetComponent(2) ;
	valZ = coord->GetComponent(3) ;
	for (i=1 ; i<=size ; i++)
	    val->At(i) = f->Evaluate(valX->At(i),valY->At(i),valZ->At(i),t,icf) ;
    }
}

