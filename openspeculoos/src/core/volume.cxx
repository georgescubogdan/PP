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

// Volume.cxx

#include "core/volume.hxx"
#include "core/face.hxx"
#include "core/edge.hxx"
#include "core/vertex.hxx"
#include "core/parent.hxx"
#include "core/point.hxx"
#include "core/distrib.hxx"
#include "core/mesh.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include "core/parallel.hxx"
#include <stdio.h>

static int ione = IONE ;


//_________________________________ Volume ____________________________________


Volume :: Volume ()
    : Element ()
{
    // Empty constructor.

    InitializeDistributions() ;
}


Volume :: ~Volume ()
{
    // Destructor.
}


Point* Volume :: GetBarycenter ()
{
    // A debugging tool.
    // Returns a new point: the barycenter of the receiver.

    return Eval(ZERO,ZERO,ZERO) ;
}


Element* Volume :: GetSubelement (int i)
{
    // Returns the 'i'-th face of the receiver.

    return GetFace(i) ;
}


void Volume :: PrintPointers ()
{
    // A debugging tool.
    // Prints the receiver (pointers only) on standard output.

    int i ;

    printf("Volume %x with faces:\n",this) ;
    for (i=1 ; i<=GetNbFaces() ; i++)
	GetFace(i)->PrintPointers() ;
}


//__________________________________ Hexa _____________________________________


int Hexa :: numberOfFaces = 6 ;             // a static attribute


Hexa :: Hexa (Quad* le, Quad* ri, Quad* fr, Quad* re, Quad* bo, Quad* to)
    : Volume ()
{
    // Constructor. Initializes the receiver to a hexaheron with the 6 faces 
    // 'le', 'ri', 'fr', 're', 'bo' and 'to'.
    //
    // These faces must be given in the following order: left, right, front, 
    // rear, bottom, top. 
    // Each of these faces can be oriented in any way.

    ref(le) ;
    ref(ri) ;
    ref(fr) ;
    ref(re) ;
    ref(bo) ;
    ref(to) ;

    faces = new Vector<Quad*>(numberOfFaces) ;
    faces->At(LEFT)   = le ;
    faces->At(RIGHT)  = ri ;
    faces->At(FRONT)  = fr ;
    faces->At(REAR)   = re ;
    faces->At(BOTTOM) = bo ;
    faces->At(TOP)    = to ;

    rotations         = NULL ;
    yaws              = NULL ;

    CheckConsistentFaces() ;
}


Hexa :: ~Hexa ()
{
    // Destructor.

    int i ;

    for (i=1 ; i<=faces->GetSize() ; i++)
	unref(faces->At(i)) ;
    delete faces ;

    delete rotations ;
    delete yaws ;
}


void Hexa :: CheckConsistentFaces ()
{
    // Issues an error message if the faces of the receiver are not properly 
    // connected.
    // Includes determining their orientation (rotation and direction).

    Face *left, *right, *front, *rear, *bottom, *top ;
    Edge *south, *north, *west, *east ;

    delete rotations ;                      // usually is already NULL
    delete yaws ;                           // usually is already NULL
    rotations = new Vector<int>(numberOfFaces) ;
    yaws      = new Vector<int>(numberOfFaces) ;
    rotations->SetValues(0) ;
    yaws->SetValues(0) ;

    left   = faces->At(LEFT) ;
    right  = faces->At(RIGHT) ;
    front  = faces->At(FRONT) ;
    rear   = faces->At(REAR) ;
    bottom = faces->At(BOTTOM) ;
    top    = faces->At(TOP) ;
  

    // 1) LEFT FACE

    south = left->GetEdge(SOUTH) ;           // edges of left face
    north = left->GetEdge(NORTH) ;
    west  = left->GetEdge(WEST) ;
    east  = left->GetEdge(EAST) ;

    if (bottom->Has(south)) {                // rotation = 0 degree
	rotations->At(LEFT) = 0 ;
	if (rear->Has(west)) {
	    yaws->At(LEFT) = 1 ;
	    if (! front->Has(east))
		Wrong(1,"left") ;
	}
	else if (front->Has(west)) {
	    yaws->At(LEFT) = 2 ;
	    if (! rear->Has(east))
		Wrong(2,"left") ;
	}
	else
	    Wrong(3,"left") ;
	if (! top->Has(north))
	    Wrong(4,"left") ;
    }

    else if (front->Has(south)) {            // rotation = 90 degrees
	rotations->At(LEFT) = 1 ;
	if (bottom->Has(west)) {
	    yaws->At(LEFT) = 1 ;
	    if (! top->Has(east))
		Wrong(5,"left") ;
	}
	else if (top->Has(west)) {
	    yaws->At(LEFT) = 2 ;
	    if (! bottom->Has(east))
		Wrong(6,"left") ;
	}
	else
	    Wrong(7,"left") ;
	if (! rear->Has(north))
	    Wrong(8,"left") ;
    }

    else if (top->Has(south)) {              // rotation = 180 degrees
	rotations->At(LEFT) = 2 ;
	if (front->Has(west)) {
	    yaws->At(LEFT) = 1 ;
	    if (! rear->Has(east))
		Wrong(9,"left") ;
	}
	else if (rear->Has(west)) {
	    yaws->At(LEFT) = 2 ;
	    if (! front->Has(east))
		Wrong(10,"left") ;
	}
	else
	    Wrong(11,"left") ;
	if (! bottom->Has(north))
	    Wrong(12,"left") ;
    }

    else if (rear->Has(south)) {             // rotation = 270 degrees
	rotations->At(LEFT) = 3 ;
	if (top->Has(west)) {
	    yaws->At(LEFT) = 1 ;
	    if (! bottom->Has(east))
		Wrong(13,"left") ;
	}
	else if (bottom->Has(west)) {
	    yaws->At(LEFT) = 2 ;
	    if (! top->Has(east))
		Wrong(14,"left") ;
	}
	else
	    Wrong(15,"left") ;
	if (! front->Has(north))
	    Wrong(16,"left") ;
    }

    else
	Error("Hexa::HasConsistentFaces","left face not connected!") ;
    
    // 2) RIGHT FACE

    south = right->GetEdge(SOUTH) ;          // edges of right face
    north = right->GetEdge(NORTH) ;
    west  = right->GetEdge(WEST) ;
    east  = right->GetEdge(EAST) ;

    if (bottom->Has(south)) {                // rotation = 0 degree
	rotations->At(RIGHT) = 0 ;
	if (front->Has(west)) {
	    yaws->At(RIGHT) = 1 ;
	    if (! rear->Has(east))
		Wrong(1,"right") ;
	}
	else if (rear->Has(west)) {
	    yaws->At(RIGHT) = 2 ;
	    if (! front->Has(east))
		Wrong(2,"right") ;
	}
	else
	    Wrong(3,"right") ;
	if (! top->Has(north))
	    Wrong(4,"right") ;
    }

    else if (rear->Has(south)) {             // rotation = 90 degrees
	rotations->At(RIGHT) = 1 ;
	if (bottom->Has(west)) {
	    yaws->At(RIGHT) = 1 ;
	    if (! top->Has(east))
		Wrong(5,"right") ;
	}
	else if (top->Has(west)) {
	    yaws->At(RIGHT) = 2 ;
	    if (! bottom->Has(east))
		Wrong(6,"right") ;
	}
	else
	    Wrong(7,"right") ;
	if (! front->Has(north))
	    Wrong(8,"right") ;
    }

    else if (top->Has(south)) {              // rotation = 180 degrees
	rotations->At(RIGHT) = 2 ;
	if (rear->Has(west)) {
	    yaws->At(RIGHT) = 1 ;
	    if (! front->Has(east))
		Wrong(9,"right") ;
	}
	else if (front->Has(west)) {
	    yaws->At(RIGHT) = 2 ;
	    if (! rear->Has(east))
		Wrong(10,"right") ;
	}
	else
	    Wrong(11,"right") ;
	if (! bottom->Has(north))
	    Wrong(12,"right") ;
    }

    else if (front->Has(south)) {            // rotation = 270 degrees
	rotations->At(RIGHT) = 3 ;
	if (top->Has(west)) {
	    yaws->At(RIGHT) = 1 ;
	    if (! bottom->Has(east))
		Wrong(13,"right") ;
	}
	else if (bottom->Has(west)) {
	    yaws->At(RIGHT) = 2 ;
	    if (! top->Has(east))
		Wrong(14,"right") ;
	}
	else
	    Wrong(15,"right") ;
	if (! rear->Has(north))
	    Wrong(16,"right") ;
    }

    else
	Error("Hexa::HasConsistentFaces","right face not connected!") ;
    
    // 3) FRONT FACE

    south = front->GetEdge(SOUTH) ;          // edges of front face
    north = front->GetEdge(NORTH) ;
    west  = front->GetEdge(WEST) ;
    east  = front->GetEdge(EAST) ;

    if (bottom->Has(south)) {                // rotation = 0 degree
	rotations->At(FRONT) = 0 ;
	if (left->Has(west)) {
	    yaws->At(FRONT) = 1 ;
	    if (! right->Has(east))
		Wrong(1,"front") ;
	}
	else if (right->Has(west)) {
	    yaws->At(FRONT) = 2 ;
	    if (! left->Has(east))
		Wrong(2,"front") ;
	}
	else
	    Wrong(3,"front") ;
	if (! top->Has(north))
	    Wrong(4,"front") ;
    }

    else if (right->Has(south)) {            // rotation = 90 degrees
	rotations->At(FRONT) = 1 ;
	if (bottom->Has(west)) {
	    yaws->At(FRONT) = 1 ;
	    if (! top->Has(east))
		Wrong(5,"front") ;
	}
	else if (top->Has(west)) {
	    yaws->At(FRONT) = 2 ;
	    if (! bottom->Has(east))
		Wrong(6,"front") ;
	}
	else
	    Wrong(7,"front") ;
	if (! left->Has(north))
	    Wrong(8,"front") ;
    }

    else if (top->Has(south)) {              // rotation = 180 degrees
	rotations->At(FRONT) = 2 ;
	if (right->Has(west)) {
	    yaws->At(FRONT) = 1 ;
	    if (! left->Has(east))
		Wrong(9,"front") ;
	}
	else if (left->Has(west)) {
	    yaws->At(FRONT) = 2 ;
	    if (! right->Has(east))
		Wrong(10,"front") ;
	}
	else
	    Wrong(11,"front") ;
	if (! bottom->Has(north))
	    Wrong(12,"front") ;
    }

    else if (left->Has(south)) {             // rotation = 270 degrees
	rotations->At(FRONT) = 3 ;
	if (top->Has(west)) {
	    yaws->At(FRONT) = 1 ;
	    if (! bottom->Has(east))
		Wrong(13,"front") ;
	}
	else if (bottom->Has(west)) {
	    yaws->At(FRONT) = 2 ;
	    if (! top->Has(east))
		Wrong(14,"front") ;
	}
	else
	    Wrong(15,"front") ;
	if (! right->Has(north))
	    Wrong(16,"front") ;
    }

    else
	Error("Hexa::HasConsistentFaces","front face not connected!") ;
    
    // 4) REAR FACE

    south = rear->GetEdge(SOUTH) ;           // edges of rear face
    north = rear->GetEdge(NORTH) ;
    west  = rear->GetEdge(WEST) ;
    east  = rear->GetEdge(EAST) ;

    if (bottom->Has(south)) {                // rotation = 0 degree
	rotations->At(REAR) = 0 ;
	if (right->Has(west)) {
	    yaws->At(REAR) = 1 ;
	    if (! left->Has(east))
		Wrong(1,"rear") ;
	}
	else if (left->Has(west)) {
	    yaws->At(REAR) = 2 ;
	    if (! right->Has(east))
		Wrong(2,"rear") ;
	}
	else
	    Wrong(3,"rear") ;
	if (! top->Has(north))
	    Wrong(4,"rear") ;
    }

    else if (left->Has(south)) {             // rotation = 90 degrees
	rotations->At(REAR) = 1 ;
	if (bottom->Has(west)) {
	    yaws->At(REAR) = 1 ;
	    if (! top->Has(east))
		Wrong(5,"rear") ;
	}
	else if (top->Has(west)) {
	    yaws->At(REAR) = 2 ;
	    if (! bottom->Has(east))
		Wrong(6,"rear") ;
	}
	else
	    Wrong(7,"rear") ;
	if (! right->Has(north))
	    Wrong(8,"rear") ;
    }

    else if (top->Has(south)) {              // rotation = 180 degrees
	rotations->At(REAR) = 2 ;
	if (left->Has(west)) {
	    yaws->At(REAR) = 1 ;
	    if (! right->Has(east))
		Wrong(9,"rear") ;
	}
	else if (right->Has(west)) {
	    yaws->At(REAR) = 2 ;
	    if (! left->Has(east))
		Wrong(10,"rear") ;
	}
	else
	    Wrong(11,"rear") ;
	if (! bottom->Has(north))
	    Wrong(12,"rear") ;
    }

    else if (right->Has(south)) {            // rotation = 270 degrees
	rotations->At(REAR) = 3 ;
	if (top->Has(west)) {
	    yaws->At(REAR) = 1 ;
	    if (! bottom->Has(east))
		Wrong(13,"rear") ;
	}
	else if (bottom->Has(west)) {
	    yaws->At(REAR) = 2 ;
	    if (! top->Has(east))
		Wrong(14,"rear") ;
	}
	else
	    Wrong(15,"rear") ;
	if (! left->Has(north))
	    Wrong(16,"rear") ;
    }

    else
	Error("Hexa::HasConsistentFaces","rear face not connected!") ;
    
    // 5) BOTTOM FACE

    south = bottom->GetEdge(SOUTH) ;         // edges of bottom face
    north = bottom->GetEdge(NORTH) ;
    west  = bottom->GetEdge(WEST) ;
    east  = bottom->GetEdge(EAST) ;

    if (front->Has(south)) {                 // rotation = 0 degree
	rotations->At(BOTTOM) = 0 ;
	if (right->Has(west)) {
	    yaws->At(BOTTOM) = 1 ;
	    if (! left->Has(east))
		Wrong(1,"bottom") ;
	}
	else if (left->Has(west)) {
	    yaws->At(BOTTOM) = 2 ;
	    if (! right->Has(east))
		Wrong(2,"bottom") ;
	}
	else
	    Wrong(3,"bottom") ;
	if (! rear->Has(north))
	    Wrong(4,"bottom") ;
    }

    else if (left->Has(south)) {             // rotation = 90 degrees
	rotations->At(BOTTOM) = 1 ;
	if (front->Has(west)) {
	    yaws->At(BOTTOM) = 1 ;
	    if (! rear->Has(east))
		Wrong(5,"bottom") ;
	}
	else if (rear->Has(west)) {
	    yaws->At(BOTTOM) = 2 ;
	    if (! front->Has(east))
		Wrong(6,"bottom") ;
	}
	else
	    Wrong(7,"bottom") ;
	if (! right->Has(north))
	    Wrong(8,"bottom") ;
    }

    else if (rear->Has(south)) {             // rotation = 180 degrees
	rotations->At(BOTTOM) = 2 ;
	if (left->Has(west)) {
	    yaws->At(BOTTOM) = 1 ;
	    if (! right->Has(east))
		Wrong(9,"bottom") ;
	}
	else if (right->Has(west)) {
	    yaws->At(BOTTOM) = 2 ;
	    if (! left->Has(east))
		Wrong(10,"bottom") ;
	}
	else
	    Wrong(11,"bottom") ;
	if (! front->Has(north))
	    Wrong(12,"bottom") ;
    }

    else if (right->Has(south)) {            // rotation = 270 degrees
	rotations->At(BOTTOM) = 3 ;
	if (rear->Has(west)) {
	    yaws->At(BOTTOM) = 1 ;
	    if (! front->Has(east))
		Wrong(13,"bottom") ;
	}
	else if (front->Has(west)) {
	    yaws->At(BOTTOM) = 2 ;
	    if (! rear->Has(east))
		Wrong(14,"bottom") ;
	}
	else
	    Wrong(15,"bottom") ;
	if (! left->Has(north))
	    Wrong(16,"bottom") ;
    }

    else
	Error("Hexa::HasConsistentFaces","bottom face not connected!") ;
    
    // 6) TOP FACE

    south = top->GetEdge(SOUTH) ;            // edges of top face
    north = top->GetEdge(NORTH) ;
    west  = top->GetEdge(WEST) ;
    east  = top->GetEdge(EAST) ;

    if (front->Has(south)) {                 // rotation = 0 degree
	rotations->At(TOP) = 0 ;
	if (left->Has(west)) {
	    yaws->At(TOP) = 1 ;
	    if (! right->Has(east))
		Wrong(1,"top") ;
	}
	else if (right->Has(west)) {
	    yaws->At(TOP) = 2 ;
	    if (! left->Has(east))
		Wrong(2,"top") ;
	}
	else
	    Wrong(3,"top") ;
	if (! rear->Has(north))
	    Wrong(4,"top") ;
    }

    else if (right->Has(south)) {            // rotation = 90 degrees
	rotations->At(TOP) = 1 ;
	if (front->Has(west)) {
	    yaws->At(TOP) = 1 ;
	    if (! rear->Has(east))
		Wrong(5,"top") ;
	}
	else if (rear->Has(west)) {
	    yaws->At(TOP) = 2 ;
	    if (! front->Has(east))
		Wrong(6,"top") ;
	}
	else
	    Wrong(7,"top") ;
	if (! left->Has(north))
	    Wrong(8,"top") ;
    }

    else if (rear->Has(south)) {             // rotation = 180 degrees
	rotations->At(TOP) = 2 ;
	if (right->Has(west)) {
	    yaws->At(TOP) = 1 ;
	    if (! left->Has(east))
		Wrong(9,"top") ;
	}
	else if (left->Has(west)) {
	    yaws->At(TOP) = 2 ;
	    if (! right->Has(east))
		Wrong(10,"top") ;
	}
	else
	    Wrong(11,"top") ;
	if (! front->Has(north))
	    Wrong(12,"top") ;
    }

    else if (left->Has(south)) {             // rotation = 270 degrees
	rotations->At(TOP) = 3 ;
	if (rear->Has(west)) {
	    yaws->At(TOP) = 1 ;
	    if (! front->Has(east))
		Wrong(13,"top") ;
	}
	else if (front->Has(west)) {
	    yaws->At(TOP) = 2 ;
	    if (! rear->Has(east))
		Wrong(14,"top") ;
	}
	else
	    Wrong(15,"top") ;
	if (! right->Has(north))
	    Wrong(16,"top") ;
    }

    else
	Error("Hexa::HasConsistentFaces","top face not connected!") ;
}

void Hexa :: CopyAddValuesFromVolume (Volume* elemX, 
				      int indx, int indy, int icx, int icy,
				      CopyAddMode cam, InterpMode im, 
				      NonConform nc, ReceiveMode rm) 
{
    // Volume from volume version of method CopyAddValuesFrom.

# ifdef REQUIRE
    Require("same object", elemX == this) ;
    Require("correct processor", Mpi_rank == procNumber) ;
    InMethod("Hexa::CopyAddValuesFromVolume") ;
# endif
 
    ElementaryField *eX, *eY, *eWork ;

    eY = GetElementaryField(indy) ;
    eX = elemX->GetElementaryField(indx) ;

    if (eX->HasSameInterpolationAs(eY)) {
    
	// 1. Conformal case

	if (cam == ADD_CAM)
	    eY->GetComponent(icy)->Add(eX->GetComponent(icx)) ;
	else
	    eY->GetComponent(icy)->CopyFrom(eX->GetComponent(icx)) ;
    }

    else {

	// 2. Non-conformal case

	eWork = eY->GetWork() ;

	if (im == NORMAL_IM)
	    eWork->CopyInterpolateFrom(eX,icx,1) ;
	else
	    eWork->CopyInterpolateTFrom(eX,icx,1) ;

	if (cam == ADD_CAM)
	    eY->GetComponent(icy)->Add(eWork->GetComponent(1)) ;
	else
	    eY->GetComponent(icy)->CopyFrom(eWork->GetComponent(1)) ; // copies twice

	ElementaryField::Retrieve(eWork) ;
    }
}


void Hexa :: CopyAddValuesFromFace (Face* elemX,
                                    int indx, int indy, int icx, int icy,
                                    CopyAddMode cam, InterpMode im,
                                    NonConform nc, ReceiveMode rm)
{
    // Volume-from-face version of method CopyAddValuesFrom.
    // No interpolation, extreme values are chosen. This means implicitly that
    // the collocation scheme of the 'indy'-th elementary field of the receiver
    // is defined until the border.

# ifdef REQUIRE
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    Require("no transposition if mortar", nc==MORTAR_NC || im!=TRANSPOSED_IM) ;
    InMethod("Hexa::CopyAddValuesFromEdge(elemX,indx,indy,icx,icy,cam,im,nc") ;
# endif

    ElementaryField     *eX, *eY, *eWork ;
    int                 dirKsi, dirEta, pY ;
    //  CopyAddMode         foo = cam ;       // just to avoid a compilation warning!

    // 1. Determine which 2 directions out of 3 of y are relevant

    GetDirectionsOf(elemX,dirKsi,dirEta) ;

    // 2. Copy values; do interpolation if required

    eY = GetElementaryField(indy) ;
    eX = elemX->GetElementaryField(indx) ;

    pY = procNumber ;

    if (eX->GetParentElement(1) == eY->GetParentElement(dirKsi) &&
	eX->GetParentElement(2) == eY->GetParentElement(dirEta)) {

	// 2a. Conformal case

	Inject(elemX,eX,eY,icx,icy,cam,rm) ;

	if (Mpi_rank == pY && rm == SIZE_RM)    // case handled in the call to
	    return ;                              // Inject herabove
    }

    else {

	// 2b. Nonconformal case: reinterpolate eX (2D) into eWork (2D); eWork has
	//     same interpolation as eY in the 2 involved directions

	if (procNumber != elemX->GetProcessorNumber())
	    NotImplemented("Hexa::CopyAddValuesFromFace",1) ;

	eWork = eX->GetWork() ;
	eWork->SetParentElement(eY->GetParentElement(dirKsi),1) ;
	eWork->SetParentElement(eY->GetParentElement(dirEta),2) ;
	if (nc == MORTAR_NC)                         // border = eWork, mortar = eX
	    elemX->ApplyQ(eWork,eX,1,icx,im) ;         // eWork = Q eX
	else {
	    if (im == TRANSPOSED_IM)
		eWork->CopyInterpolateTFrom(eX,icx,1) ;
	    else
		eWork->CopyInterpolateFrom(eX,icx,1) ;
	}

	Inject(elemX,eWork,eY,1,icy,cam,rm) ;

	if (Mpi_rank == pY && rm == SIZE_RM)    // case handled in the call to
	    return ;                              // Inject herabove

	ElementaryField::Retrieve(eWork) ;
    }
}


Point* Hexa :: Eval (real r, real s, real t)
{
    // Returns a new point: the receiver evaluated at (r,s,t).

# ifdef REQUIRE
    Require("'r' is in [-1,1}", r >= -1. && r <= 1.) ;
    Require("'s' is in [-1,1}", s >= -1. && s <= 1.) ;
    Require("'t' is in [-1,1}", t >= -1. && t <= 1.) ;
    InMethod("Hexa::Eval(r,s,t)") ;
# endif

    Face  *left, *right, *front, *rear, *bottom, *top ;
    Point *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *p9, *p10, *p11, *p12, *p13,
        *p14, *p15, *p16, *p17, *p18, *p19, *p20, *p21, *p22, *p23, *p24,
        *p25, *p26 ;
    real  a, b, phi1r, phi2r, phi1s, phi2s, phi1t, phi2t ;


    phi1r  = Phi1(r) ;
    phi2r  = Phi2(r) ;
    phi1s  = Phi1(s) ;
    phi2s  = Phi2(s) ;
    phi1t  = Phi1(t) ;
    phi2t  = Phi2(t) ;
    left   = GetFace(LEFT) ;
    right  = GetFace(RIGHT) ;
    front  = GetFace(FRONT) ;
    rear   = GetFace(REAR) ;
    bottom = GetFace(BOTTOM) ;
    top    = GetFace(TOP) ;

    // contributions of the 6 faces

    VolumeToFace(M1,s,t,LEFT,a,b) ;
    p1 = left->Eval(a,b)->Multiply(phi1r) ; 
    VolumeToFace(P1,s,t,RIGHT,a,b) ;
    p2 = right->Eval(a,b)->Multiply(phi2r) ;
    VolumeToFace(r,M1,t,FRONT,a,b) ;
    p3 = front->Eval(a,b)->Multiply(phi1s) ;
    VolumeToFace(r,P1,t,REAR,a,b) ;
    p4 = rear->Eval(a,b)->Multiply(phi2s) ;
    VolumeToFace(r,s,M1,BOTTOM,a,b) ;
    p5 = bottom->Eval(a,b)->Multiply(phi1t) ;
    VolumeToFace(r,s,P1,TOP,a,b) ;
    p6 = top->Eval(a,b)->Multiply(phi2t) ;

    // contributions of the 12 edges

    VolumeToFace(M1,M1,t,LEFT,a,b) ;
    p7  = left->Eval(a,b)->Multiply(-phi1r*phi1s) ;
    VolumeToFace(M1,P1,t,LEFT,a,b) ;
    p8  = left->Eval(a,b)->Multiply(-phi1r*phi2s) ;
    VolumeToFace(P1,M1,t,RIGHT,a,b) ;
    p9  = right->Eval(a,b)->Multiply(-phi2r*phi1s) ;
    VolumeToFace(P1,P1,t,RIGHT,a,b) ;
    p10 = right->Eval(a,b)->Multiply(-phi2r*phi2s) ;
    VolumeToFace(r,M1,M1,FRONT,a,b) ;
    p11 = front->Eval(a,b)->Multiply(-phi1s*phi1t) ;
    VolumeToFace(r,M1,P1,FRONT,a,b) ;
    p12 = front->Eval(a,b)->Multiply(-phi1s*phi2t) ;
    VolumeToFace(r,P1,M1,REAR,a,b) ;
    p13 = rear->Eval(a,b)->Multiply(-phi2s*phi1t) ;
    VolumeToFace(r,P1,P1,REAR,a,b) ;
    p14 = rear->Eval(a,b)->Multiply(-phi2s*phi2t) ;
    VolumeToFace(M1,s,M1,BOTTOM,a,b) ;
    p15 = bottom->Eval(a,b)->Multiply(-phi1t*phi1r) ;
    VolumeToFace(P1,s,M1,BOTTOM,a,b) ;
    p16 = bottom->Eval(a,b)->Multiply(-phi1t*phi2r) ;
    VolumeToFace(M1,s,P1,TOP,a,b) ;
    p17 = top->Eval(a,b)->Multiply(-phi2t*phi1r) ;
    VolumeToFace(P1,s,P1,TOP,a,b) ;
    p18 = top->Eval(a,b)->Multiply(-phi2t*phi2r) ;

    // contributions of the 8 vertices

    VolumeToFace(M1,M1,M1,LEFT,a,b) ;
    p19 = left->Eval(a,b)->Multiply(phi1r*phi1s*phi1t) ;
    VolumeToFace(M1,M1,P1,LEFT,a,b) ;
    p20 = left->Eval(a,b)->Multiply(phi1r*phi1s*phi2t) ;
    VolumeToFace(M1,P1,M1,LEFT,a,b) ;
    p21 = left->Eval(a,b)->Multiply(phi1r*phi2s*phi1t) ;
    VolumeToFace(M1,P1,P1,LEFT,a,b) ;
    p22 = left->Eval(a,b)->Multiply(phi1r*phi2s*phi2t) ;
    VolumeToFace(P1,M1,M1,RIGHT,a,b) ;
    p23 = right->Eval(a,b)->Multiply(phi2r*phi1s*phi1t) ;
    VolumeToFace(P1,M1,P1,RIGHT,a,b) ;
    p24 = right->Eval(a,b)->Multiply(phi2r*phi1s*phi2t) ;
    VolumeToFace(P1,P1,M1,RIGHT,a,b) ;
    p25 = right->Eval(a,b)->Multiply(phi2r*phi2s*phi1t) ;
    VolumeToFace(P1,P1,P1,RIGHT,a,b) ;
    p26 = right->Eval(a,b)->Multiply(phi2r*phi2s*phi2t) ;



    p1->Add(p2)->Add(p3)->Add(p4)->Add(p5)->Add(p6)->Add(p7)->Add(p8)->Add(p9)
	->Add(p10)->Add(p11)->Add(p12)->Add(p13)->Add(p14)->Add(p15)->Add(p16)
	->Add(p17)->Add(p18)->Add(p19)->Add(p20)->Add(p21)->Add(p22)->Add(p23)
	->Add(p24)->Add(p25)->Add(p26) ;

    unref(p2)  ; unref(p3)  ; unref(p4)  ; unref(p5)  ; unref(p6)  ;
    unref(p7)  ; unref(p8)  ; unref(p9)  ; unref(p10) ; unref(p11) ;
    unref(p12) ; unref(p13) ; unref(p14) ; unref(p15) ; unref(p16) ;
    unref(p17) ; unref(p18) ; unref(p19) ; unref(p20) ; unref(p21) ;
    unref(p22) ; unref(p23) ; unref(p24) ; unref(p25) ; unref(p26) ;

    return p1 ;
}


ElementaryField* Hexa :: Extract (int indy, int icy, Face* face,
                                  ReceiveMode rm) 
{
    // Returns a 1-component elementary field on 'face': the restriction to 
    // 'face' of the 'icy'-th component of the 'indy'-th elementary field of the
    // receiver.
    //
    // In order to avoid dynamic allocations, the answer is a work elementary
    // field, so get rid of it by means of method 'Retrieve' rather than the
    // 'delete' operator.

# ifdef REQUIRE
    Require("face belongs to volume", Has(face)) ;
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(face->GetProcessorNumber())) ;
    InMethod("Volume::Extract(indy,icy,face,rm)") ;
# endif

    ElementaryField     *eY, *answer ;
    ParentElement       *parent1, *parent2 ;
    RealVector          *valY, *valA ;
    int                 dir1, dir2, k, location, np1, np2, np3, startY,
	startA, size, jump, pX, pY ;

    // 1. Determine the position of the face on the volume

    location = SearchLocationOf(face) ;

    // 2. Create the returned field (2D) with the adequate parents.

    pY  = procNumber ;
    pX  = face->GetProcessorNumber() ;
    eY  = GetElementaryField(indy) ;
    np1 = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    np2 = eY->GetParentElement(2)->GetNbCollocationPoints() ;
    np3 = eY->GetParentElement(3)->GetNbCollocationPoints() ;

    if (Mpi_rank == pX && rm == SIZE_RM) {
	if (location == LEFT || location == RIGHT)
	    size = np2 * np3 ;
	else if (location == FRONT || location == REAR)
	    size = np1 * np3 ;
	else  // location == BOTTOM || location == TOP) 
	    size = np1 * np2 ;
	Mpi_receiveBufferSizes[pY] += size ;
	return NULL ;
    }

    GetDirectionsOfFace(location,dir1,dir2) ;
    parent1 = eY->GetParentElement(dir1) ;
    parent2 = eY->GetParentElement(dir2) ;
    answer  = eY->GetWork(parent1,parent2) ;                       // 2D

    // 3. Extract the values of the face into 'answer'. The values are not stored
    //    in the correct order, due to the 8 possible orientations of the face.
    //    For example, if 'face' is the front or the rear face, the values are 
    //    ordered according to the (r,t) axes of the receiver, not according to
    //    the axes of 'face'.

    if (Mpi_rank == pY)
	valY = eY->GetComponent(icy) ;                               // 3D
    if (Mpi_rank == pX)
	valA = answer->GetComponent(1) ;                             // 2D

    if (location == LEFT || location == RIGHT) {
	if (location == LEFT)
	    startY = 1 ;
	else
	    startY = np1 ;
	startA = 1 ;
	jump   = np1*np2 ;
	for (k=1 ; k<=np3 ; k++) {
	    if (pX == pY)                                     // same processor
		DCOPY(&np2, &valY->At(startY), &np1, &valA->At(startA), &ione) ;
	    else {                                            // different processors
		if (Mpi_rank == pY)                             // for 3D processor
		    Mpi_SendToBuffer(np2, &valY->At(startY), np1, pX) ;
		else                                            // for 2D processor
		    Mpi_ReceiveFromBuffer(np2, &valA->At(startA), ione, pY) ;
	    }
	    startY += jump ;
	    startA += np2 ;
	}
    }

    else if (location == FRONT || location == REAR) {
	if (location == FRONT)
	    startY = 1 ;
	else
	    startY = (np2-1)*np1 + 1 ;
	startA = 1 ;
	jump   = np1*np2 ;
	for (k=1 ; k<=np3 ; k++) {
	    if (pX == pY)                                     // same processor
		DCOPY(&np1, &valY->At(startY), &ione, &valA->At(startA), &ione) ;
	    else {                                            // different processors
		if (Mpi_rank == pY)                             // for 3D processor
		    Mpi_SendToBuffer(np1, &valY->At(startY), ione, pX) ;
		else                                            // for 2D processor
		    Mpi_ReceiveFromBuffer(np1, &valA->At(startA), ione, pY) ;
	    }
	    startY += jump ;
	    startA += np1 ;
	}
    }

    else { // location == BOTTOM || location == TOP) 
	if (location == BOTTOM)
	    startY = 1 ;
	else
	    startY = (np3-1)*np1*np2 + 1 ;
	size = np1*np2 ;
	if (pX == pY)                                       // same processor
	    DCOPY(&size, &valY->At(startY), &ione, &valA->At(1), &ione) ;
	else {                                              // different processors
	    if (Mpi_rank == pY)                               // for 3D processor
		Mpi_SendToBuffer(size, &valY->At(startY), ione, pX) ;
	    else                                              // for 2D processor
		Mpi_ReceiveFromBuffer(size, &valA->At(1), ione, pY) ;
	}
    }

    // 4. Reorder the values of 'answer' by taking orientation into account

    if (Mpi_rank == pX)
	ReindexToFace(face,answer,1) ;

    return answer ;
}


Mesh3D* Hexa :: GenerateMesh () 
{
    // Returns a new mesh obtained by refining the receiver using the receiver's
    // distributions.
    // In order account for the fact faces and edges may not be straight/plane,
    // Gordon-Hall transfinite interpolations are used to generate vertices,
    // edges and faces. Exception: interior edges are generated as straight 
    // lines. Exterior edges are generated by transfinite interpolations (ie,
    // they are curved), in order to represent curved boundaries accurately.

    Vector<Quad*>   *facesRS, *facesST, *facesTR ;
    Vector<Edge*>   *edgesR, *edgesS, *edgesT ;
    Vector<Vertex*> *vertices ;
    Mesh3D          *mesh ;
    Hexa            *volume ;
    Quad            *face, *left, *right, *front, *rear, *bottom, *top ;
    Edge            *edge, *south, *north, *west, *east ;
    Vertex          *v1, *v2 ;
    Distribution    *distribR, *distribS, *distribT ;
    Point           *newPoint ;
    real            r, r1, r2, s, s1, s2, t, t1, t2, a1, a2, b1, b2 ;
    int             i, j, k, ii, jj, kk, pos, index, npR, npS, npT, 
                    nsegR, nsegS, nsegT, faceNum ;

    // Generation of the vertices

    Printf("\nStart mesh generation ...") ;

    distribR = GetDistribution(1) ;
    distribS = GetDistribution(2) ;
    distribT = GetDistribution(3) ;
    npR      = distribR->GetSize() ;
    npS      = distribS->GetSize() ;
    npT      = distribT->GetSize() ;
    vertices = new Vector<Vertex*>(npR*npS*npT) ;
    index    = 1 ;
    for (k=1 ; k<=npT ; k++) {
	t = distribT->At(k) ;
	for (j=1 ; j<=npS ; j++) {
	    s = distribS->At(j) ;
	    for (i=1 ; i<=npR ; i++) {
		r                   = distribR->At(i) ;
		newPoint            = Eval(r,s,t) ;
		vertices->At(index) = new Vertex(newPoint) ;
		index++ ;
	    }
	}
    }

    // Generation of the edges
    // a1,a2,b1,b2 = 2 of the 3 coordinates (r,s,t), taking face orientation into
    // account 

    nsegR  = npR - 1 ;
    nsegS  = npS - 1 ;
    nsegT  = npT - 1 ;
    edgesR = new Vector<Edge*>(nsegR*npS*npT) ;
    edgesS = new Vector<Edge*>(nsegS*npT*npR) ;
    edgesT = new Vector<Edge*>(nsegT*npR*npS) ;

    index = 1 ;
    for (k=1 ; k<=npT ; k++) {                      // edges in direction r
	t = distribT->At(k) ;
	for (j=1 ; j<=npS ; j++) {
	    s = distribS->At(j) ;
	    for (i=1 ; i<=nsegR ; i++) {
		pos = (k-1)*npR*npS + (j-1)*npR + i ;
		v1  = vertices->At(pos) ;
		v2  = vertices->At(pos+1) ;
		r1  = distribR->At(i) ;
		r2  = distribR->At(i+1) ;
		if (k==1 || k==npT || j==1 || j==npS) {   // external edge
		    if (k == 1)
			faceNum = BOTTOM ;
		    else if (k == npT) 
			faceNum = TOP ;
		    else if (j == 1)
			faceNum = FRONT ;
		    else if (j == npS)
			faceNum = REAR ;
		    VolumeToFace(r1,s,t,faceNum,a1,b1) ;
		    VolumeToFace(r2,s,t,faceNum,a2,b2) ;
		    edge = faces->At(faceNum)->ExtractEdgeBetween(min(a1,a2),
								  max(a1,a2),min(b1,b2),max(b1,b2)) ;
		    edge->SetVerticesTo(v1,v2) ;
		}
		else                                      // internal edge
		    edge = new Line(v1,v2) ;
		edgesR->At(index++) = edge ;
	    }
	}
    }

    index = 1 ;
    for (i=1 ; i<=npR ; i++) {                      // edges in direction s
	r = distribR->At(i) ;
	for (k=1 ; k<=npT ; k++) {
	    t = distribT->At(k) ;
	    for (j=1 ; j<=nsegS ; j++) {
		pos = (k-1)*npR*npS + (j-1)*npR + i ;
		v1  = vertices->At(pos) ;
		v2  = vertices->At(pos+npR) ;
		s1  = distribS->At(j) ;
		s2  = distribS->At(j+1) ;
		if (i==1 || i==npR || k==1 || k==npT) {   // external edge
		    if (i == 1)
			faceNum = LEFT ;
		    else if (i == npR)
			faceNum = RIGHT ;
		    else if (k == 1)
			faceNum = BOTTOM ;
		    else if (k == npT) 
			faceNum = TOP ;
		    VolumeToFace(r,s1,t,faceNum,a1,b1) ;
		    VolumeToFace(r,s2,t,faceNum,a2,b2) ;
		    edge = faces->At(faceNum)->ExtractEdgeBetween(min(a1,a2),
								  max(a1,a2),min(b1,b2),max(b1,b2)) ;
		    edge->SetVerticesTo(v1,v2) ;
		}
		else                                      // internal edge
		    edge = new Line(v1,v2) ;
		edgesS->At(index++) = edge ;
	    }
	}
    }

    index = 1 ;
    for (j=1 ; j<=npS ; j++) {                      // edges in direction t
	s = distribS->At(j) ;
	for (i=1 ; i<=npR ; i++) {
	    r = distribR->At(i) ;
	    for (k=1 ; k<=nsegT ; k++) {
		pos = (k-1)*npR*npS + (j-1)*npR + i ;
		v1  = vertices->At(pos) ;
		v2  = vertices->At(pos+npR*npS) ;
		t1  = distribT->At(k) ;
		t2  = distribT->At(k+1) ;
		if (j==1 || j==npS || i==1 || i==npR) {   // external edge
		    if (j == 1)
			faceNum = FRONT ;
		    else if (j == npS)
			faceNum = REAR ;
		    else if (i == 1)
			faceNum = LEFT ;
		    else if (i == npR) 
			faceNum = RIGHT ;
		    VolumeToFace(r,s,t1,faceNum,a1,b1) ;
		    VolumeToFace(r,s,t2,faceNum,a2,b2) ;
		    edge = faces->At(faceNum)->ExtractEdgeBetween(min(a1,a2),
								  max(a1,a2),min(b1,b2),max(b1,b2)) ;
		    edge->SetVerticesTo(v1,v2) ;
		}
		else                                      // internal edge
		    edge = new Line(v1,v2) ;
		edgesT->At(index++) = edge ;
	    }
	}
    }

    // Generation of the faces
 
    facesRS = new Vector<Quad*>(nsegR*nsegS*npT) ;
    facesST = new Vector<Quad*>(nsegS*nsegT*npR) ;
    facesTR = new Vector<Quad*>(nsegT*nsegR*npS) ;

    index = 1 ;
    for (k=1 ; k<=npT ; k++) {                      // faces in the r-s planes
	ii = (k-1)*nsegR*npS + 1 ;
	for (j=1 ; j<=nsegS ; j++) {
	    jj = (k-1)*nsegS + j ;
	    for (i=1 ; i<=nsegR ; i++) {
		south = edgesR->At(ii) ;
		north = edgesR->At(ii+nsegR) ;
		west  = edgesS->At(jj) ;
		east  = edgesS->At(jj+nsegS*npT) ;
		face  = new Quad(south,north,west,east) ;
		facesRS->At(index++) = face ;
		ii++ ;
		jj += nsegS*npT ;
	    }
	}
    }

    index = 1 ;
    for (i=1 ; i<=npR ; i++) {                      // faces in the s-t planes
	jj = (i-1)*nsegS*npT + 1 ;
	for (k=1 ; k<=nsegT ; k++) {
	    kk = (i-1)*nsegT + k ;
	    for (j=1 ; j<=nsegS ; j++) {
		south = edgesS->At(jj) ;
		north = edgesS->At(jj+nsegS) ;
		west  = edgesT->At(kk) ;
		east  = edgesT->At(kk+nsegT*npR) ;
		face  = new Quad(south,north,west,east) ;
		facesST->At(index++) = face ;
		jj++ ;
		kk += nsegT*npR ;
	    }
	}
    }

    index = 1 ;
    for (j=1 ; j<=npS ; j++) {                      // faces in the t-r planes
	kk = (j-1)*nsegT*npR + 1 ;
	for (i=1 ; i<=nsegR ; i++) {
	    ii = (j-1)*nsegR + i ;
	    for (k=1 ; k<=nsegT ; k++) {
		south = edgesT->At(kk) ;
		north = edgesT->At(kk+nsegT) ;
		west  = edgesR->At(ii) ;
		east  = edgesR->At(ii+nsegR*npS) ;
		face  = new Quad(south,north,west,east) ;
		facesTR->At(index++) = face ;
		kk++ ;
		ii += nsegR*npS ;
	    }
	}
    }

    // Generation of the volumes

    mesh = new Mesh3D() ;
    for (k=1 ; k<=nsegT ; k++)                      // volumes
	for (j=1 ; j<=nsegS ; j++)
	    for (i=1 ; i<=nsegR ; i++) {
		ii = (i-1)*nsegS*nsegT + (k-1)*nsegS + j ;
		jj = (j-1)*nsegT*nsegR + (i-1)*nsegT + k ;
		kk = (k-1)*nsegR*nsegS + (j-1)*nsegR + i ;
		left   = facesST->At(ii) ;
		right  = facesST->At(ii+nsegS*nsegT) ;
		front  = facesTR->At(jj) ;
		rear   = facesTR->At(jj+nsegT*nsegR) ;
		bottom = facesRS->At(kk) ;
		top    = facesRS->At(kk+nsegR*nsegS) ;
		volume = new Hexa(left,right,front,rear,bottom,top) ;
		mesh->PutWithoutCheck(volume) ;
		unref(volume) ;
	    }

    // cleaning
    for (i=1 ; i<=vertices->GetSize() ; i++)
	unref(vertices->At(i)) ;
    delete vertices ;

    for (i=1 ; i<=edgesR->GetSize() ; i++)
	unref(edgesR->At(i)) ;
    delete edgesR ;

    for (i=1 ; i<=edgesS->GetSize() ; i++)
	unref(edgesS->At(i)) ;
    delete edgesS ;

    for (i=1 ; i<=edgesT->GetSize() ; i++)
	unref(edgesT->At(i)) ;
    delete edgesT ;

    for (i=1 ; i<=facesRS->GetSize() ; i++)
	unref(facesRS->At(i)) ;
    delete facesRS ;

    for (i=1 ; i<=facesST->GetSize() ; i++)
	unref(facesST->At(i)) ;
    delete facesST ;

    for (i=1 ; i<=facesTR->GetSize() ; i++)
	unref(facesTR->At(i)) ;
    delete facesTR ;

    Printf(" done\n") ;

    return mesh ;
}


void Hexa :: GetDirectionsOf (Face* f, int& dir1, int& dir2)
{
    // Initializes 'dir1' and 'dir2' to two values equal to 1,2,3, the directions
    // of 'f' in the receiver.

# ifdef REQUIRE
    Require("face 'f' belongs to hexahedron", Has(f)) ;
    InMethod("Hexa::GetDirectionsOf(i,dir1,dir2)") ;
# endif

    GetDirectionsOfFace(faces->GetIndexOf((Quad*)f),dir1,dir2) ;
}


void Hexa :: GetDirectionsOfFace (int i, int& dir1, int& dir2)
{
    // Initializes 'dir1' and 'dir2' to two values in {1,2,3}, the directions
    // of the 'i'-th face of the receiver.

# ifdef REQUIRE
    Require("valid index 'i'", i >= 1 && i <= numberOfFaces) ;
    InMethod("Hexa::GetDirectionsOfFace(i,dir1,dir2)") ;
# endif

    int rot ;

    rot = rotations->At(i) ;

    if (i == FRONT || i == REAR) {
	if (rot == 0 || rot == 2) {
	    dir1 = 1 ;
	    dir2 = 3 ;
	}
	else {
	    dir1 = 3 ;
	    dir2 = 1 ;
	}
    }

    else if (i == LEFT || i == RIGHT) {
	if (rot == 0 || rot == 2) {
	    dir1 = 2 ;
	    dir2 = 3 ;
	}
	else {
	    dir1 = 3 ;
	    dir2 = 2 ;
	}
    }

    else if (i == BOTTOM || i == TOP) {
	if (rot == 0 || rot == 2) {
	    dir1 = 1 ;
	    dir2 = 2 ;
	}
	else {
	    dir1 = 2 ;
	    dir2 = 1 ;
	}
    }
}
  

Face* Hexa :: GetFace (int i)
{
    // Returns the i-th face of the receiver.
    // Remark: C++ is a shit: the return value should be allowed to be declared
    //         as a Quad, not just as a Face.

# ifdef REQUIRE
    Require("valid face number 'i'", i >= 1 && i <= GetNbFaces()) ;
    InMethod("Hexa::GetFace(i)") ;
# endif

    return faces->At(i) ;
}


void Hexa :: GetSetInteriorDof (int index, int iop, int idof, real& val)
{
    // Gets ('iop'=1) or sets ('iop'=0) the 'idof'-th interior degree of freedom 
    // to the 'index'-th elementary field of the receiver from/to 'val'.
    // An interior degree of freedom is a degree of freedom that is not involved
    // in a continuity constraint. For Gauss-Lobatto collocation grids, these are
    // precisely the points located at the interior of the element. For Gauss
    // grids, there should not not be any distinction between interior or usual
    // degree of freedom. However, it is assumed that discretization is based on
    // Gauss-Lobatto.

# ifdef REQUIRE
    Require("'index' is a valid field number", HasElementaryField(index)) ;
    Require("valid parent element 1",  
	    GetElementaryField(index)->GetParentElement(1)->IsDefinedOnBorders()) ;
    Require("valid parent element 2",  
	    GetElementaryField(index)->GetParentElement(2)->IsDefinedOnBorders()) ;
    Require("valid parent element 3",  
	    GetElementaryField(index)->GetParentElement(3)->IsDefinedOnBorders()) ;
    InMethod("Hexa::GetSetInteriorDof(index,iop,idof,val)") ;
# endif

    ElementaryField *eY ;
    RealVector      *valY ;
    int             i, j, k, ic, np1, np2, npi1, npi2, npi3, loc, size, 
	layerSize ;

    eY   = GetElementaryField(index) ;
    np1  = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    np2  = eY->GetParentElement(2)->GetNbCollocationPoints() ;
    npi1 = eY->GetParentElement(1)->GetNbInteriorCollocationPoints() ;
    npi2 = eY->GetParentElement(2)->GetNbInteriorCollocationPoints() ;
    npi3 = eY->GetParentElement(3)->GetNbInteriorCollocationPoints() ;

    size      = npi1 * npi2 * npi3 ;
    ic        = (idof-1)/size + 1 ;
    idof      = (idof-1)%size + 1 ;
    layerSize = npi1 * npi2 ;
    k         = (idof-1)/layerSize + 1 ;
    idof      = (idof-1)%layerSize + 1 ;
    j         = (idof-1)/npi1 + 1 ;
    i         = (idof-1)%npi1 + 1 ;
    loc       = k*np1*np2 + j*np1 + i+1 ;

    valY = eY->GetComponent(ic) ;
    if (iop)
	valY->At(loc) = val ;
    else
	val = valY->At(loc) ;
}


real Hexa :: GetTypicalLength (int dir)
{
    // Returns a rough approximation of the mean length of the receiver in the
    // local 'dir'-th direction.

# ifdef REQUIRE
    Require("valid argument 'dir'", dir == 1 || dir == 2 || dir == 3) ;
    InMethod("Hexa::GetTypicalLength(dir)") ;
# endif

    Point *p1, *p2 ;
    real  answer ;

    if (dir == 1) {
	p1 = Eval(-ONE,ZERO,ZERO) ;
	p2 = Eval( ONE,ZERO,ZERO) ;
    }
    else if (dir == 2) {
	p1 = Eval(ZERO,-ONE,ZERO) ;
	p2 = Eval(ZERO, ONE,ZERO) ;
    }
    else {
	p1 = Eval(ZERO,ZERO,-ONE) ;
	p2 = Eval(ZERO,ZERO, ONE) ;
    }
  
    p2->Subtract(p1) ;
    answer = p2->GetNorm() ;

    delete p1 ;
    delete p2 ;

    return answer ;
}


boolean Hexa :: Has (Element* elem)
{
    // Returns True if 'elem' is one of the faces or one of the edges or one of
    // the vertices (same pointer) of the receiver.

    Face *face ;
    int  i, j, n, m ;

    if (elem->IsFace())
	return faces->Has((Quad*)elem) ;

    else if (elem->IsEdge()) {
	n = GetNbFaces() ;
	for (i=1 ; i<=n ; i++)
	    if (faces->At(i)->Has((Edge*)elem))
		return true ;
	return false ;
    }

    else if (elem->IsVertex()) {
	n = GetNbFaces() ;
	for (i=1 ; i<=n ; i++) {
	    face = faces->At(i) ;
	    m    = face->GetNbEdges() ;
	    for (j=1 ; j<=m ; j++)
		if (face->GetEdge(j)->Has((Vertex*)elem))
		    return true ;
	    return false ;
	}
    }

    return false ;
}


void Hexa :: Inject (Face* elemX, ElementaryField* eX, ElementaryField* eY,
                     int icx, int icy, CopyAddMode cam, ReceiveMode rm)
{
    // Copies or adds the 'icx'-th component of 'eX' to the 'icy'-th component of
    // 'eY'.
    // 'eX' belongs to the face 'elemX' of the receiver.
    // 'eY' belongs to the receiver.
    // 'eX' is assumed to have the same interpolation as 'eY'. 
    // Accounts for the orientation of 'elemX' in the receiver.

# ifdef REQUIRE
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Hexa::Inject(elemX,eX,eY,icx,icy") ;
# endif

    ElementaryField *eWork ;
    RealVector      *valY, *valW ;
    real            one=ONE ;
    int             location, startY, startW, jump, np1, np2, np3, k, ic, size,
	pX, pY, mode ;

    pY = procNumber ;
    pX = elemX->GetProcessorNumber() ; 

    // 1. Reindex the values in eX (2D).
    //    The values of eX cannot be simply injected into eY, due to the 8 
    //    possible orientations of the face 'elemX'.
    //    For example, if 'face' is the front or the rear face, the values are 
    //    ordered according to the (ksi,eta) axes of the face, and thus must be
    //    reordered according to the (r,t) axes of the receiver.
    //    Avoids the use of a work field (2D) if the 2 systems of axes coincide.

    eWork = NULL ;
    if (Mpi_rank == pX) {
	eWork = ReindexToVolume(elemX,eX,icx) ;
	if (eWork == NULL) {          // same axes: no creation of eWork, so use eX
	    eWork = eX ;
	    ic    = icx ;
	}
	else                          // eWork is a 1-component field
	    ic = 1 ;
    }

    // 2. Inject the values of eWork (2D) into eY (3D)

    if (Mpi_rank == pY)
	valY = eY->GetComponent(icy) ;                                  // 3D
    if (Mpi_rank == pX)
	valW = eWork->GetComponent(ic) ;                                // 2D

    np1 = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    np2 = eY->GetParentElement(2)->GetNbCollocationPoints() ;
    np3 = eY->GetParentElement(3)->GetNbCollocationPoints() ;

    location = SearchLocationOf(elemX) ;

    if (Mpi_rank == pY && rm == SIZE_RM) {
	if (location == LEFT || location == RIGHT)
	    size = np2 * np3 ;
	else if (location == FRONT || location == REAR)
	    size = np1 * np3 ;
	else  // location == BOTTOM || location == TOP 
	    size = np1 * np2 ;
	Mpi_receiveBufferSizes[pX] += size ;
	return ;
    }

    if (cam == COPY_CAM)
	mode = 0 ;
    else
	mode = 1 ;

    if (location == LEFT || location == RIGHT) {
	if (location == LEFT)
	    startY = 1 ;
	else
	    startY = np1 ;
	startW = 1 ;
	jump   = np1*np2 ;
	for (k=1 ; k<=np3 ; k++) {
	    if (pX == pY) {
		if (cam == COPY_CAM)
		    DCOPY(&np2, &valW->At(startW), &ione, &valY->At(startY), &np1) ;
		else
		    DAXPY(&np2, &one, &valW->At(startW), &ione, &valY->At(startY), &np1);
	    }
	    else {
		if (Mpi_rank == pX)
		    Mpi_SendToBuffer(np2, &valW->At(startW), ione, pY) ;
		else 
		    Mpi_ReceiveFromBuffer(np2, &valY->At(startY), np1, pX, mode) ;
	    }
	    startY += jump ;
	    startW += np2 ;
	}
    }

    else if (location == FRONT || location == REAR) {
	if (location == FRONT)
	    startY = 1 ;
	else
	    startY = (np2-1)*np1 + 1 ;
	startW = 1 ;
	jump   = np1*np2 ;
	for (k=1 ; k<=np3 ; k++) {
	    if (pX == pY) {
		if (cam == COPY_CAM)
		    DCOPY(&np1, &valW->At(startW), &ione, &valY->At(startY), &ione) ;
		else
		    DAXPY(&np1, &one, &valW->At(startW), &ione, &valY->At(startY),&ione);
	    }
	    else {
		if (Mpi_rank == pX)
		    Mpi_SendToBuffer(np1, &valW->At(startW), ione, pY) ;
		else 
		    Mpi_ReceiveFromBuffer(np1, &valY->At(startY), ione, pX, mode) ;
	    }
	    startY += jump ;
	    startW += np1 ;
	}
    }

    else { // location == BOTTOM || location == TOP) 
	if (location == BOTTOM)
	    startY = 1 ;
	else
	    startY = (np3-1)*np1*np2 + 1 ;
	size = np1*np2 ;
	if (pX == pY) {
	    if (cam == COPY_CAM)
		DCOPY(&size, &valW->At(1), &ione, &valY->At(startY), &ione) ;
	    else
		DAXPY(&size, &one, &valW->At(1), &ione, &valY->At(startY), &ione) ;
	}
	else {
	    if (Mpi_rank == pX)
		Mpi_SendToBuffer(size, &valW->At(1), ione, pY) ;
	    else 
		Mpi_ReceiveFromBuffer(size, &valY->At(startY), ione, pX, mode) ;
	}
    }

    if (Mpi_rank == pX && eWork != eX)
	ElementaryField::Retrieve(eWork) ;
}


void Hexa :: Print ()
{
    // Prints the receiver on standard output.

    int i ;

    printf ("Hexa:\n") ;

    if (rotations)
	printf ("  rotations:  %d %d %d %d %d %d\n",
		rotations->At(1),rotations->At(2),rotations->At(3),
		rotations->At(4),rotations->At(5),rotations->At(6)) ;
    else
	printf ("  rotations:  NULL\n") ;

    if (yaws)
	printf ("  yaws: %d %d %d %d %d %d\n\n",
		yaws->At(1),yaws->At(2),yaws->At(3),
		yaws->At(4),yaws->At(5),yaws->At(6)) ;
    else
	printf ("  yaws: NULL\n\n") ;

    for (i=1 ; i<=6 ; i++) {
	printf ("Face %d: ",i) ; 
	faces->At(i)->Print() ;
	printf ("\n") ;
    }
}


void Hexa :: ReindexToFace (Face* elemX, ElementaryField* eX, int icx)
{
    // The values of the 'icx'-th component of the elementary field 'eX' of the
    // face 'elemX' of the receiver are assumed to be ordered according to the
    // axes of the receiver. 
    // Reorders these values according to the axes of the face 'elemX', not to 
    // the axes of the receiver.

# ifdef REQUIRE
    Require("face belongs to volume", Has(elemX)) ;
    Require("valid index 'icx'", icx >= 1 && icx <= eX->GetNbComponents()) ;
    Require("correct processor", Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Hexa::ReindexToFace(elemX,eX,icx)") ;
# endif

    ParentElement *parent1, *parent2 ;
    RealVector    *valX, *valT ;
    real          sOld, sNew, tOld, tNew ;
    int           i, location, np1, np2, size, startT, startX, incrX ;

    // Find the orientation of the face on the receiver

    location = SearchLocationOf(elemX) ;
    sOld     = THIRD ;                          // any value
    tOld     = HALF ;                           // any value different from sOld
    VolumeFaceToFace(sOld,tOld,location,sNew,tNew) ;

    // Avoid copies if volume and face have coincident axes

    if (sNew == sOld && tNew == tOld) 
	// case 1: special case: nothing to do
	return ;

    // Copy the values into a temporary storage

    valX = eX->GetComponent(icx) ;
    valT = new RealVector(valX->GetSize()) ;
    valT->CopyFrom(valX) ;

    // Reorder the values (see hand-written document by YDP)

    parent1 = eX->GetParentElement(1) ;
    parent2 = eX->GetParentElement(2) ;
    np1     = parent1->GetNbCollocationPoints() ;
    np2     = parent2->GetNbCollocationPoints() ;
    size    = np1 * np2 ;
    startT  = 1 ;

    if (tNew == tOld) {                                 // and thus sNew = -sOld
	// case 2: first row becomes first row (reversed), etc
	startX =  1 ;
	incrX  = -1 ;
	for (i=1 ; i<=np2 ; i++) {
	    DCOPY(&np1, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
	    startT += np1 ;
	    startX += np1 ;
	}
    }

    else if (tNew == -sOld) {
	if (sNew == tOld) {
	    // case 3: first row becomes first column (reversed), etc
	    startX = 1 ;
	    incrX  = -np1 ;
	    for (i=1 ; i<=np1 ; i++) {
		DCOPY(&np2, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np2 ;
		startX++ ;
	    }
	}
	else {                                            // sNew = -tOld
	    // case 4: first row becomes last column (reversed), etc
	    startX =  np1 ; 
	    incrX  = -np1 ;
	    for (i=1 ; i<=np1 ; i++) {
		DCOPY(&np2, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np2 ;
		startX-- ;
	    }
	}
    }

    else if (tNew == -tOld) {
	if (sNew == -sOld) {
	    // case 5: first row becomes last row (reversed), etc
	    startX = size - np1 + 1 ;
	    incrX  = -1 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np1 ;
		startX -= np1 ;
	    }
	}
	else {                                            // sNew = sOld
	    // case 6: first row becomes last row, etc
	    startX = size - np1 + 1 ;
	    incrX  = 1 ; 
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np1 ;
		startX -= np1 ;
	    }
	}
    }

    else if (tNew == sOld) {     
	if (sNew == -tOld) {
	    // case 7: first row becomes last column, etc
	    startX = np1 ;
	    incrX  = np1 ;
	    for (i=1 ; i<=np1 ; i++) {
		DCOPY(&np2, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np2 ;
		startX-- ;
	    }
	}
	else {                                            // sNew = tOld
	    // case 8: first row becomes first column, etc
	    startX = 1 ;
	    incrX  = np1 ;
	    for (i=1 ; i<=np1 ; i++) {
		DCOPY(&np2, &valT->At(startT), &ione, &valX->At(startX), &incrX) ;
		startT += np2 ;
		startX++ ;
	    }
	}
    }

    else
	Error("Volume::ReindexToFace(elemX,eX,icx)","wrong value of (sNew,tNew)") ;

    delete valT ;
}


ElementaryField* Hexa :: ReindexToVolume (Face* elemX, ElementaryField* eX, 
                                          int icx)
{
    // Returns a new elementary field with 1 component, equal to the 'icx'-th
    // component of the elementary field 'eX' supposed to belong to the face
    // 'elemX' of the receiver; however the new field has values ordered
    // according to the axes of the receiver, rather than according to the axes
    // of the face 'elemX'.
    // Returns NULL if the axes of 'elemX' coincide with those of the receiver.

# ifdef REQUIRE
    Require("face belongs to volume", Has(elemX)) ;
    Require("valid index 'icx'", icx >= 1 && icx <= eX->GetNbComponents()) ;
    Require("correct processor", Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Hexa::ReindexToVolume(elemX,eX,icx)") ;
# endif

    ElementaryField *answer ;
    ParentElement   *parent1, *parent2 ;
    RealVector      *valX, *valA ;
    real            sOld, sNew, tOld, tNew ;
    int             i, location, np1, np2, size, startX, startA, incrA ;

    // Find the orientation of the face on the receiver

    location = SearchLocationOf(elemX) ;
    sOld     = THIRD ;                           // any value
    tOld     = HALF ;                            // any value different from sOld
    VolumeFaceToFace(sOld,tOld,location,sNew,tNew) ;

    // Avoid copies if volume and face have coincident axes

    if (sNew == sOld && tNew == tOld) 
	// case 1: special case: nothing to do
	return NULL ;

    // Reorder the values (see hand-written document by YDP)

    parent1 = eX->GetParentElement(1) ;
    parent2 = eX->GetParentElement(2) ;
    np1     = parent1->GetNbCollocationPoints() ;
    np2     = parent2->GetNbCollocationPoints() ;
    valX    = eX->GetComponent(icx) ;
    size    = np1 * np2 ;
    startX  = 1 ;

    if (tNew == tOld) {                                 // and thus sNew = -sOld
	// case 2: first row (reversed) becomes first row, etc
	answer = ElementaryField::GetWork(parent1,parent2) ;
	valA   = answer->GetComponent(1) ;
	startA =  1 ;
	incrA  = -1 ;
	for (i=1 ; i<=np2 ; i++) {
	    DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
	    startX += np1 ;
	    startA += np1 ;
	}
    }

    else if (tNew == -sOld) {
	answer = ElementaryField::GetWork(parent2,parent1) ;
	valA   = answer->GetComponent(1) ;
	if (sNew == tOld) {
	    // case 3: first row becomes last column, etc
	    //         unlike method ReindexToFace!
	    startA = np2 ;
	    incrA  = np2 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA++ ;
	    }
	}
	else {                                            // sNew = -tOld
	    // case 4: first row becomes last column (reversed), etc
	    startA =  np2 ; 
	    incrA  = -np2 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA-- ;
	    }
	}
    } 

    else if (tNew == -tOld) {
	answer = ElementaryField::GetWork(parent1,parent2) ;
	valA   = answer->GetComponent(1) ;
	if (sNew == -sOld) {
	    // case 5: first row becomes last row (reversed), etc
	    startA = size - np1 + 1 ;
	    incrA  = -1 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA -= np1 ;
	    }
	}
	else {                                            // sNew = sOld
	    // case 6: first row becomes last row, etc
	    startA = size - np1 + 1 ;
	    incrA  = 1 ; 
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA -= np1 ;
	    }
	}
    }

    else {                                              // tNew = sOld 
	answer = ElementaryField::GetWork(parent2,parent1) ;
	valA   = answer->GetComponent(1) ;
	if (sNew == -tOld) {
	    // case 7: first row becomes first column (reversed), etc
	    startA =  1 ;
	    incrA  = -np2 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA-- ;
	    }
	}
	else {                                            // sNew = tOld
	    // case 8: first row becomes first column, etc
	    startA = 1 ;
	    incrA  = np2 ;
	    for (i=1 ; i<=np2 ; i++) {
		DCOPY(&np1, &valX->At(startX), &ione, &valA->At(startA), &incrA) ;
		startX += np1 ;
		startA++ ;
	    }
	}
    }

    return answer ;
}


int Hexa :: SearchLocationOf (Element* face)
{
    // Returns the position index of 'face' in the receiver. Returns 0 if
    // 'face' does not belong to the receiver.

# ifdef REQUIRE
    Require("'face' is a face", face->IsFace()) ;
    Require("'face' is a quad", ((Face*)face)->IsQuad()) ;
    InMethod("Hexa::SearchLocationOf(face)") ;
# endif

    return faces->GetIndexOf((Quad*)face) ;
}


void Hexa :: SetToCoordinates (int index)
{
    // Sets the elementary field 'index' to the coordinates of the receiver.
    // The receiver initializes as many components as it finds.
    // In every component the coordinates are stored as follows: first, the
    // Bottom layer t=-1, next the next layer, etc, until the Top layer t=1.
    // At every layer the coordinates are stored row by row: first the row s=-1,
    // next the next row, etc, until the last row s=1. On every row the 
    // coordinates are stored from r=-1 to r=1 (the r and s local coordinates
    // denote here those of the volume, not of a face!).

# ifdef REQUIRE
    Require("number of components <= 3",
	    GetElementaryField(index)->GetNbComponents() <= 3) ;
    InMethod("Face::SetToCoordinates") ;
# endif

    ElementaryField *eY ;
    Point           *point ;
    RealVector      *collocationPoints1, *collocationPoints2, 
	*collocationPoints3 ;
    real            r, s, t ;
    int             i, j, k, icomp, np1, np2, np3, ncomp, pos ;

    if (Mpi_rank != procNumber)
	return ;

    eY                 = GetElementaryField(index) ;
    ncomp              = eY->GetNbComponents() ;
    collocationPoints1 = eY->GetParentElement(1)->GetCollocationPoints();
    collocationPoints2 = eY->GetParentElement(2)->GetCollocationPoints();
    collocationPoints3 = eY->GetParentElement(3)->GetCollocationPoints();
    np1                = collocationPoints1->GetSize() ;
    np2                = collocationPoints2->GetSize() ;
    np3                = collocationPoints3->GetSize() ;
    pos                = 0 ;
    for (k=1 ; k<=np3 ; k++) {
	t = collocationPoints3->At(k) ;
	for (j=1 ; j<=np2 ; j++) {
	    s = collocationPoints2->At(j) ;
	    for (i=1 ; i<=np1 ; i++) {
		r     = collocationPoints1->At(i) ;
		point = Eval(r,s,t) ;
		pos++ ;
		for (icomp=1; icomp<=ncomp ; icomp++)
		    eY->GetComponent(icomp)->At(pos) = point->GetCoordinate(icomp) ;
		unref(point) ;
	    }
	}
    }
}


void Hexa :: Substitute (Face* source, Face* target)
{
    // Assuming that 'source' is one of the faces (same object) of the receiver,
    // replaces 'source' by 'target' in the face list of the receiver.
    // 'source' and 'target' are expected to be instances of Quad, not just of 
    // Face; they are expected to be geometrically coincident; they need not be
    // oriented identically.

# ifdef REQUIRE
    Require("'source' is a Quad", source->IsQuad()) ;
    Require("'target' is a Quad", target->IsQuad()) ;
    Require("'source' and 'target' are geometrically coincident",
            target->IsIdenticalTo(source)) ;
    Require("'source' belongs to hexahedron", 
            faces->GetIndexOf((Quad*)source) != 0) ;
    InMethod("Hexa::Substitute(source,target)") ;
# endif

    Quad *src, *trg ;
    Edge *south, *west ;
    int  i, index ;

    ref(target) ;

    src = (Quad*)source ;
    trg = (Quad*)target ;

    // substitute the face ('source' with 'target')
    index = faces->GetIndexOf(src) ;
    faces->At(index) = trg ;

    // update the rotation of the face
    south = target->GetEdge(SOUTH) ;
    for (i=1 ; i<=4 ; i++)
	if (source->GetEdge(i)->IsIdenticalTo(south))
	    break ; 
    rotations->At(index) += i-1 ;
    if (rotations->At(index) >= 4)
	rotations->At(index) -= 4 ;              // e.g., if rot=6, then rot=2

    // update the direction of the face
    west = target->GetEdge(WEST) ;
    i   += 2 ;                                 // position of western face of src
    i    = i % 4 ;                             // e.g., if i=6, then i=2
    if (! west->IsIdenticalTo(src->GetEdge(i))) {  // if trg and src have opposi-
	yaws->At(index) += 1 ;                   //                  te directions
	if (yaws->At(index) == 3)
	    yaws->At(index) = 1 ;
    }

    // update the other faces (they must share edges with 'target', no longer
    // with 'source')
    for (i=1 ; i<=numberOfFaces ; i++)
	if (i != index)
	    GetFace(i)->SubstituteCommonEdges(target) ;

    unref(source) ;
}


void Hexa :: VolumeFaceToFace (real A, real B, int f, real& a, real& b)
{
    // Converts abscissae from volume's face viewpoint to face viewpoint.
    // (A,B) contains the abscissae of a point on the 'f'-th face of the 
    // receiver, in the axes of the receiver. For example, if f=REAR, then (A,B)
    // are the abscissae (r,t).
    // Stores in (a,b) the equivalent of (A,B) in the orientation of the face.

# ifdef REQUIRE
    Require("'A' is in [-1,1]", A >= -1. && A <= 1.) ;
    Require("'B' is in [-1,1]", B >= -1. && B <= 1.) ;
    Require("valid argument 'f'", f >= 1 && f <= numberOfFaces) ;
    InMethod("Hexa::VolumeFaceToFace(A,B,f,a,b)") ;
# endif

    real r, s, t ;

    if (f == LEFT) {
	r = ONE_MINUS ;
	s = A ;
	t = B ;
    }

    else if (f == RIGHT) {
	r = ONE ;
	s = A ;
	t = B ;
    }

    else if (f == FRONT) {
	r = A ;
	s = ONE_MINUS ;
	t = B ;
    }

    else if (f == REAR) {
	r = A ;
	s = ONE ;
	t = B ;
    }

    else if (f == BOTTOM) {
	r = A ;
	s = B ;
	t = ONE_MINUS ;
    }

    else {// f == TOP
	r = A ;
	s = B ;
	t = ONE ;
    }

    VolumeToFace(r,s,t,f,a,b) ;
}


void Hexa :: VolumeToFace (real r, real s, real t, int f, real& a, real& b)
{
    // Converts abscissae from volume viewpoint to face viewpoint.
    // (r,s,t) contains the abscissae of a point on the surface; 'f' is the 
    // number (betwen 1 and 6) of the face on which that point stands; the
    // corresponding blocked abscissa (r, or s, or t) must be equal to -1 or 1.
    // Stores in 'a' and 'b' the equivalent 2 free abscissae (as opposed to the
    // third, blocked, abscissa), accounting for the orientation of the face in
    // the receiver.
    // Example: if eg f=FRONT (s must be equal to -1.); a and b are such that 
    // the following two messages are equivalent:
    //    . point = this->Eval(r,s,t) ;
    //    . point = faces->At(f)->Eval(a,b) ;

# ifdef REQUIRE
    Require("'r' is in [-1,1]", r >= -1. && r <= 1.) ;
    Require("'s' is in [-1,1]", s >= -1. && s <= 1.) ;
    Require("'t' is in [-1,1]", t >= -1. && t <= 1.) ;
    Require("valid argument 'f'", f >= 1 && f <= numberOfFaces) ;
    Require("on left face, 'r' == -1.",  f != LEFT   || r == -1.) ;
    Require("on right face, 'r' == 1.",  f != RIGHT  || r ==  1.) ;
    Require("on front face, 's' == -1.", f != FRONT  || s == -1.) ;
    Require("on rear face, 's' == 1.",   f != REAR   || s ==  1.) ;
    Require("on bottom face, 't' == -1.",f != BOTTOM || t == -1.) ;
    Require("on top face, 't' == 1.",    f != TOP    || t ==  1.) ;
    InMethod("Hexa::VolumeToFace(r,s,t,f,a,b)") ;
# endif

    int rot, yaw ;

    rot = rotations->At(f) ;
    yaw = yaws->At(f) ;
    a   = 100. ;   // dummy values, for debugging
    b   = 100. ;

    if (f == LEFT) {
	if (rot==0 && yaw==1) {
	    a = -s ;
	    b = +t ;
	}
	else if (rot==0 && yaw==2) {
	    a = +s ;
	    b = +t ;
	}
	else if (rot==1 && yaw==1) {
	    a = +t ;
	    b = +s ;
	}
	else if (rot==1 && yaw==2) {
	    a = -t ;
	    b = +s ;
	}
	else if (rot==2 && yaw==1) {
	    a = +s ;
	    b = -t ;
	}
	else if (rot==2 && yaw==2) {
	    a = -s ;
	    b = -t ;
	}
	else if (rot==3 && yaw==1) {
	    a = -t ;
	    b = -s ;
	}
	else if (rot==3 && yaw==2) {
	    a = +t ;
	    b = -s ;
	}
    }

    else if (f == RIGHT) {
	if (rot==0 && yaw==1) {
	    a = +s ;
	    b = +t ;
	}
	else if (rot==0 && yaw==2) {
	    a = -s ;
	    b = +t ;
	}
	else if (rot==1 && yaw==1) {
	    a = +t ;
	    b = -s ;
	}
	else if (rot==1 && yaw==2) {
	    a = -t ;
	    b = -s ;
	}
	else if (rot==2 && yaw==1) {
	    a = -s ;
	    b = -t ;
	}
	else if (rot==2 && yaw==2) {
	    a = +s ;
	    b = -t ;
	}
	else if (rot==3 && yaw==1) {
	    a = -t ;
	    b = +s ;
	}
	else if (rot==3 && yaw==2) {
	    a = +t ;
	    b = +s ;
	}
    }

    else if (f == FRONT) {
	if (rot==0 && yaw==1) {
	    a = +r ;
	    b = +t ;
	}
	else if (rot==0 && yaw==2) {
	    a = -r ;
	    b = +t ;
	}
	else if (rot==1 && yaw==1) {
	    a = +t ;
	    b = -r ;
	}
	else if (rot==1 && yaw==2) {
	    a = -t ;
	    b = -r ;
	}
	else if (rot==2 && yaw==1) {
	    a = -r ;
	    b = -t ;
	}
	else if (rot==2 && yaw==2) {
	    a = +r ;
	    b = -t ;
	}
	else if (rot==3 && yaw==1) {
	    a = -t ;
	    b = +r ;
	}
	else if (rot==3 && yaw==2) {
	    a = +t ;
	    b = +r ;
	}
    }

    else if (f == REAR) {
	if (rot==0 && yaw==1) {
	    a = -r ;
	    b = +t ;
	}
	else if (rot==0 && yaw==2) {
	    a = +r ;
	    b = +t ;
	}
	else if (rot==1 && yaw==1) {
	    a = +t ;
	    b = +r ;
	}
	else if (rot==1 && yaw==2) {
	    a = -t ;
	    b = +r ;
	}
	else if (rot==2 && yaw==1) {
	    a = +r ;
	    b = -t ;
	}
	else if (rot==2 && yaw==2) {
	    a = -r ;
	    b = -t ;
	}
	else if (rot==3 && yaw==1) {
	    a = -t ;
	    b = -r ;
	}
	else if (rot==3 && yaw==2) {
	    a = +t ;
	    b = -r ;
	}
    }

    else if (f == BOTTOM) {
	if (rot==0 && yaw==1) {
	    a = -r ;
	    b = +s ;
	}
	else if (rot==0 && yaw==2) {
	    a = +r ;
	    b = +s ;
	}
	else if (rot==1 && yaw==1) {
	    a = +s ;
	    b = +r ;
	}
	else if (rot==1 && yaw==2) {
	    a = -s ;
	    b = +r ;
	}
	else if (rot==2 && yaw==1) {
	    a = +r ;
	    b = -s ;
	}
	else if (rot==2 && yaw==2) {
	    a = -r ;
	    b = -s ;
	}
	else if (rot==3 && yaw==1) {
	    a = -s ;
	    b = -r ;
	}
	else if (rot==3 && yaw==2) {
	    a = +s ;
	    b = -r ;
	}
    }

    else if (f == TOP) {
	if (rot==0 && yaw==1) {
	    a = +r ;
	    b = +s ;
	}
	else if (rot==0 && yaw==2) {
	    a = -r ;
	    b = +s ;
	}
	else if (rot==1 && yaw==1) {
	    a = +s ;
	    b = -r ;
	}
	else if (rot==1 && yaw==2) {
	    a = -s ;
	    b = -r ;
	}
	else if (rot==2 && yaw==1) {
	    a = -r ;
	    b = -s ;
	}
	else if (rot==2 && yaw==2) {
	    a = +r ;
	    b = -s ;
	}
	else if (rot==3 && yaw==1) {
	    a = -s ;
	    b = +r ;
	}
	else if (rot==3 && yaw==2) {
	    a = +s ;
	    b = +r ;
	}
    }
}


void Hexa :: Wrong (int i, string face)
{
    // Issues an error message for method "HasConsistentFaces".

    printf("\nError in Hexa::HasConsistentFaces :\n") ;
    printf("       error %d when checking %s face\n\n",i,face.c_str()) ;
    Error("Hexa::HasConsistentFaces","") ;
}


Vertex* Hexa :: GetVertex (int i)
{
    if (i <= 4)
	return faces->At(5)->GetVertex(i);  // bottom face
    else
	return faces->At(6)->GetVertex(i-4);  // top face
}
