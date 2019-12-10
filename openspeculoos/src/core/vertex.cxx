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

// vertex.cxx

#include "core/vertex.hxx"
#include "core/edge.hxx"
#include "core/face.hxx"
#include "core/point.hxx"
#include "core/parent.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <stdio.h>


Vertex :: Vertex (Point* p)
    : Element ()
{
    // Constructor. Initializes the receiver to a vertex at 'p'.

    ref(p) ;
    point = p ;
    InitializeDistributions() ;
}


Vertex :: ~Vertex ()
{
    // Destructor.

    unref(point) ;
}


void Vertex :: CopyAddValuesFromVertex (Vertex* elemX, 
                                        int indx, int indy, int icx, int icy,
                                        CopyAddMode cam, InterpMode, 
                                        NonConform, ReceiveMode) 
{
    // Vertex-from-vertex version of method CopyAddValuesFrom.

# ifdef REQUIRE
    Require("same object", elemX == this) ;
    Require("correct processors", Mpi_rank == procNumber) ;
    InMethod("Vertex::CopyAddValuesFromVertex");
# endif

    RealVector *valY, *valX ;

    valY = GetElementaryField(indy)->GetComponent(icy) ;
    valX = elemX->GetElementaryField(indx)->GetComponent(icx) ;
  
    if (cam == ADD_CAM) 
	valY->At(1) += valX->At(1) ;
    else
	valY->At(1)  = valX->At(1) ;
}


void Vertex :: CopyAddValuesFromEdge (Edge* elemX, 
                                      int indx, int indy, int icx, int icy,
                                      CopyAddMode cam, InterpMode, 
                                      NonConform, ReceiveMode rm) 
{
    // Vertex-from-edge version of method CopyAddValuesFrom.

# ifdef REQUIRE
    Require("vertex belongs to edge", elemX->Has(this)) ;
    Require("parent element goes until border",
	    elemX->GetElementaryField(indx)->GetParentElement(1)->IsDefinedOnBorders());
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Vertex::CopyAddValuesFromEdge") ;
# endif

    ElementaryField *eY, *eX ;
    real            val[1] ;
    int             ival, location, pX, pY ;

    pY = procNumber ;
    pX = elemX->GetProcessorNumber() ;

    if (Mpi_rank == pY && rm == SIZE_RM) {
	Mpi_receiveBufferSizes[pX] += 1 ;
	return ;
    }

    eY       = GetElementaryField(indy) ;
    eX       = elemX->GetElementaryField(indx) ;
    location = elemX->SearchLocationOf(this) ;
    if (location == 1)                                           // left vertex
	ival = 1 ;
    else                                                         // right vertex
	ival = eX->GetNbDofs() ;

    if (pX == pY) {                                // same processor
	if (cam == ADD_CAM) 
	    eY->GetComponent(icy)->At(1) += eX->GetComponent(icx)->At(ival) ;
	else
	    eY->GetComponent(icy)->At(1)  = eX->GetComponent(icx)->At(ival) ;
    }
    else {                                         // different processors
	if (pX == Mpi_rank)                          // run by elemX's processor
	    Mpi_SendToBuffer(1, &eX->GetComponent(icx)->At(ival), 1, pY) ;
	else {                                       // run by elemY's processor
	    Mpi_ReceiveFromBuffer(1, &val[0], 1, pX) ;
	    if (cam == ADD_CAM)
		eY->GetComponent(icy)->At(1) += val[0] ;
	    else
		eY->GetComponent(icy)->At(1)  = val[0] ;
	}
    }
}


void Vertex :: CopyAddValuesFromFace (Face* elemX, 
                                      int indx, int indy, int icx, int icy,
                                      CopyAddMode cam, InterpMode, 
                                      NonConform, ReceiveMode rm) 
{
    // Vertex-from-face version of method CopyAddValuesFrom.

# ifdef REQUIRE
    Require("vertex belongs to face", elemX->Has(this)) ;
    Require("parent element 1 goes until border",
	    elemX->GetElementaryField(indx)->GetParentElement(1)->IsDefinedOnBorders());
    Require("parent element 2 goes until border",
	    elemX->GetElementaryField(indx)->GetParentElement(2)->IsDefinedOnBorders());
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Vertex::CopyAddValuesFromFace") ;
# endif

    ElementaryField *eX ;
    RealVector      *valX, *valY ;
    real            value ;
    int             i, n, np1, np2, ival, pX, pY ;
    boolean         found = false ;

    // 0. 

    pX = elemX->GetProcessorNumber() ;
    pY = procNumber ;

    if (Mpi_rank == pY && rm == SIZE_RM) {
	Mpi_receiveBufferSizes[pX] += 1 ;
	return ;
    }

    // 1. Find the position of the vertex in the face

    n = elemX->GetNbVertices() ;
    for (i=1 ; i<=n ; i++)
	if (this == elemX->GetVertex(i)) {
	    found = true ;
	    break ;
	}
    if (! found)
	Error("Vertex::CopyAddValuesFromFace","vertex not found in the face") ;

    eX  = elemX->GetElementaryField(indx) ;
    np1 = eX->GetParentElement(1)->GetNbCollocationPoints() ;
    np2 = eX->GetParentElement(2)->GetNbCollocationPoints() ;

    if (i == 1)                    // South West
	ival = 1 ; 
    else if (i == 2)               // South East
	ival = np1 ;
    else if (i == 3)               // North East
	ival = np1*np2 ;
    else                           // North West
	ival = (np1-1)*np2 + 1 ;

    // 2. Get the value from the face

    if (Mpi_rank == pX)
	valX = eX->GetComponent(icx) ;
    if (Mpi_rank == pY)
	valY = GetElementaryField(indy)->GetComponent(icy) ;

    if (pX == pY)                                        // same processor
	value = valX->At(ival) ;
    else {                                               // different processors
	if (Mpi_rank == pX)                                // for elemX's processor
	    Mpi_SendToBuffer(1, &valX->At(ival), 1, pY) ;
	else                                               // for elemY's processor
	    Mpi_ReceiveFromBuffer(1, &value, 1, pX) ;
    }
  
    // 3. Copy the value to the vertex

    if (Mpi_rank == pY) {
	if (cam == ADD_CAM) 
	    valY->At(1) += value ;
	else
	    valY->At(1)  = value ;
    }
}


Point* Vertex :: GetBarycenter ()
{
    // A debugging tool.
    // Returns a copy of the receiver's point.

    return point->Duplicate() ;
}


boolean Vertex :: IsIdenticalTo (Element* elem)
{
    // Returns True if the receiver and 'elem' coincide, else returns False.

    Vertex *v ;

    if (! elem->IsVertex())
	return false ;
    else
	v = (Vertex*) elem ;

    return (this==v || point->IsIdenticalTo(v->point)) ;
}


boolean Vertex :: IsVertex ()
{
    // Returns True, because the receiver is a vertex.

    return true ;
}


void Vertex :: Print ()
{
    // Prints the receiver on the standard output.
 
    printf ("Vertex at ") ;
    point->Print() ;
}


int Vertex :: SearchLocationOf (Element*)
{
    // Issues an error message.

    Error("Vertex::SearchLocationOf(elem)","nonsense") ;
    return 0 ;
}


void Vertex :: SetToCoordinates (int index) 
{
    // Sets the elementary field 'index' to the coordinates of the receiver.
    // The receiver initializes as many components as it finds.

    ElementaryField *eY ;
    int             i, ncomp ;

# ifdef REQUIRE
    Require("number of components <= 3", 
	    GetElementaryField(index)->GetNbComponents() <= 3) ;
    InMethod("Vertex::SetToCoordinates(index)") ;
# endif

    if (Mpi_rank != procNumber)
	return ;

    eY    = GetElementaryField(index) ;
    ncomp = eY->GetNbComponents() ;

    for (i=1; i<=ncomp ; i++) 
	eY->GetComponent(i)->At(1) = point->GetCoordinate(i) ;
}


