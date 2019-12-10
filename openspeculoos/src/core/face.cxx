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

// face.cxx

#include "core/face.hxx"


static real one  = ONE ;
static int  ione = IONE ;


//___________________________________ Face ____________________________________


Face :: Face ()
    : Element () 
{
    // Empty constructor.

    edges        = NULL ;
    orientations = NULL ;
    InitializeDistributions() ;
}


Face :: ~Face ()
    // Destructor.
{
    int i ;

    for (i=1 ; i<=edges->GetSize() ; i++)
	unref(edges->At(i)) ;
    delete edges ;

    delete orientations ;
}


Point* Face :: GetBarycenter ()
{
    // A debugging tool.
    // Returns a new point: the barycenter of the receiver.

    return Eval(ZERO,ZERO) ;
}


Edge* Face :: GetEdge (int i)
{
    // Returns the i-th edge of the receiver.

# ifdef REQUIRE
    Require("valid edge number 'i'", i >= 1 && i <= GetNbEdges()) ;
    InMethod("Face::GetEdge(i)") ;
# endif

    return edges->At(i) ;
}


Element* Face :: GetSubelement (int i)
{
    // Returns the 'i'-th edge of the receiver.

    return GetEdge(i) ;
}


boolean Face :: Has (Element* elem)
{
    // Returns True if 'elem' is one of the edges or one of the vertices (same
    // pointer) of the receiver.

    int i, n ;

    if (elem->IsEdge())
	return edges->Has((Edge*)elem) ;

    else if (elem->IsVertex()) {
	n = GetNbEdges() ;
	for (i=1 ; i<=n ; i++)
	    if (edges->At(i)->Has((Vertex*)elem))
		return true ;
	return false ;
    }
    else
	return false ;
} 


boolean Face :: IsFace ()
{
    // Returns True because the receiver is a face.

    return true ;
}


boolean Face :: IsQuad ()
{
    // Returns False because the receiver is not a Quad.

    return false ;
}


void Face :: PrintPointers ()
{
    // A debugging tool.
    // Prints on standard output the edges of the receiver (pointers only).

    int i ;

    printf("  Face %x with edges:\n",this) ;
    for (i=1 ; i<=GetNbEdges() ; i++)
	GetEdge(i)->PrintPointers() ;
}


int Face :: SearchLocationOf (Element* edge)
{
    // Returns the position index of 'edge' in the receiver. Returns 0 if
    // 'edge' does not belong to the receiver.

# ifdef REQUIRE
    Require("'edge' is an edge", edge->IsEdge()) ;
    InMethod("Face::SearchLocationOf(edge)") ;
# endif

    return edges->GetIndexOf((Edge*)edge) ;
}


void Face :: Substitute (Edge* source, Edge* target)
{
    // Substitutes edge 'source' by edge 'target'.
    // 'source' is supposed to be an edge of the receiver.
    // 'source' and 'target' are supposed to be geometrically coincident; they
    // need not point to the same vertices; they need not be oriented 
    // identically.
    // Note that several edges of the receiver are affected (must point to other
    // vertices).

# ifdef REQUIRE
    Require("source and target are identical", target->IsIdenticalTo(source)) ;
    Require("source belongs to the face", edges->GetIndexOf(source) != 0) ; 
    InMethod("Face::Substitute(source,target)") ;
# endif

    int i, j ;

    ref(target) ;

    // substitute 'source' with 'target'
    i = edges->GetIndexOf(source) ;
  
    edges->At(i) = target ;
    if (! source->HasSameOrientationAs(target)) 
	orientations->At(i) = - orientations->At(i) ;
    unref(source) ;
    // update vertices of other edges with those of 'target'
    for (j=1 ; j<=GetNbEdges() ; j++)
	if (j != i) {
	    GetEdge(j)->SubstituteCommonVertices(target) ;
	}
}


void Face :: SubstituteCommonEdges (Face* face)
{
    // For each edge e1 of the receiver that is identical to an edge e2 of 
    // 'face', makes the receiver point to e2 rather than to e1.

    Edge *e1, *e2 ;
    int  i1, i2, n1, n2 ;

    n1 = GetNbEdges () ;
    n2 = face->GetNbEdges () ;
    for (i2=1 ; i2<=n2 ; i2++) {
	e2 = face->GetEdge(i2) ;
	for (i1=1 ; i1<=n1 ; i1++) {
	    e1 = GetEdge(i1) ;
	    if (e1->IsIdenticalTo(e2))
		Substitute(e1,e2) ;
	}
    }
}


//___________________________________ Quad ____________________________________


int Quad :: numberOfEdges = 4 ; 

// initialization of a static attribute
  

boolean Quad :: curvedInteriorEdges = true ;

// initialization of a static attribute


Quad :: Quad (Edge* s, Edge* n, Edge* w, Edge* e)
    : Face ()
{
    // Constructor. Initializes the receiver to a face with the 4 edges 's', 'n',
    // 'w', 'e'.
    // Although these edges need not be oriented in any specific way, they must
    // be given in an order such that, conventionally, 's' will be termed the 
    // "South edge", 'n' the "North edge", 'w' the "West edge" and 'e' the "East
    // edge", ie, edges 1 and 2 have no common vertex (and idem for edges 3 and
    // 4). 
    // Two connected edges must point to the same vertex (not just to two 
    // geometrically coincident vertices).
    // Connectivity and convexity of the face are checked.

    ref(s) ;
    ref(n) ;
    ref(w) ;
    ref(e) ;
    edges = new Vector<Edge*>(4) ;
    edges->At(SOUTH) = s ;
    edges->At(NORTH) = n ;
    edges->At(WEST)  = w ;
    edges->At(EAST)  = e ;

    if (! HasConsistentEdges())
	Error("Quad::Quad","wrong choice of edges") ;
}


Quad :: ~Quad ()
    // Destructor.
{
}


void Quad :: ApplyQ (ElementaryField* border, ElementaryField* mortar, 
                     int icb, int icm, InterpMode im)                   
{
    // Applies the mortar operator Q:
    // - if 'im' = NORMAL_IM,     sets u(border) to  Q.u(mortar)
    // - if 'im' = TRANSPOSED_IM, sets u(mortar) to tQ.u(border)
    // Only the 'icb'-th component of 'border' and the 'icm'-th component of
    // 'mortar' are involved.
    // 'border' and 'mortar' are defined on the receiver; one or both of their
    // parent elements are different.

# ifdef REQUIRE
    Require("same nb components", border->GetNbComponents() 
	    == mortar->GetNbComponents()) ;
    Require("valid argument 'icb'", icb>=1 && icb<=border->GetNbComponents()) ;
    Require("valid argument 'icm'", icm>=1 && icm<=mortar->GetNbComponents()) ;
    InMethod("Quad::ApplyQ(border,mortar,icb,icm,im)") ;
# endif

    ElementaryField *temp ;
    ParentElement   *borderKsi, *borderEta, *mortarKsi, *mortarEta,
	*quadratureKsi, *quadratureEta ;
    int             degreeBorderKsi, degreeBorderEta,
	degreeMortarKsi, degreeMortarEta ;

    // 1. Determine the richest quadrature rules, to ensure exact integration

    borderKsi = border->GetParentElement(1) ;
    borderEta = border->GetParentElement(2) ;
    mortarKsi = mortar->GetParentElement(1) ;
    mortarEta = mortar->GetParentElement(2) ;

    degreeBorderKsi = borderKsi->GetPolynomialDegree() ;
    degreeBorderEta = borderEta->GetPolynomialDegree() ;
    degreeMortarKsi = mortarKsi->GetPolynomialDegree() ;
    degreeMortarEta = mortarEta->GetPolynomialDegree() ;

    if (degreeMortarKsi >= degreeBorderKsi)
	quadratureKsi = mortarKsi ;
    else
	quadratureKsi = borderKsi ;

    if (degreeMortarEta >= degreeBorderEta)
	quadratureEta = mortarEta ;
    else
	quadratureEta = borderEta ;

    // 2. Project on 'temp' the values of the source field ('mortar', in normal
    //    mode)
    //    Several cases are not implemented (e.g. mortar richer than border in 
    //    one direction, but poorer in the other direction). Implementing them 
    //    would require implementing methods in class TensorMatrix for mixed
    //    interpolation/transposed interpolation.

    if (im == NORMAL_IM) {                            // u(border) = Q u(mortar)
	if (quadratureKsi == mortarKsi && quadratureEta == mortarEta) {
	    temp = mortar->DuplicateEmpty(1) ;
	    temp->MultiplyByWeights(1) ;  
	    border->CopyInterpolateTFrom(temp,1,icb) ;
	    border->DivideByWeights(icb) ;
	    delete temp ;
	}
	else if (quadratureKsi == borderKsi && quadratureEta == borderEta)
	    border->CopyInterpolateFrom(mortar,icm,icb) ;
	else
	    Error("Quad::ApplyQ","case not implemented") ;
    }
    
    else {                                          // u(mortar) = tQ u(border)
	if (quadratureKsi == mortarKsi && quadratureEta == mortarEta) {
	    temp = border->DuplicateEmpty(icb) ;
	    temp->DivideByWeights(1) ;
	    mortar->CopyFrom(temp,1,icm) ;
	    mortar->MultiplyByWeights(icm) ;
	    delete temp ;
	}
	else if (quadratureKsi == borderKsi && quadratureEta == borderEta)
	    mortar->CopyInterpolateTFrom(border,icb,icm) ;
	else
	    Error("Quad::ApplyQ","case not implemented") ;
    }
}


void Quad :: CopyAddValuesFromEdge (Edge* elemX, 
                                    int indx, int indy, int icx, int icy,
                                    CopyAddMode cam, InterpMode im, 
                                    NonConform nc, ReceiveMode rm) 
{
    // Face from edge version of method CopyAddValuesFrom.
    // No interpolation, extreme values are chosen. This means implicitly that
    // the collocation scheme of the 'indy'-th elementary field of the receiver 
    // is defined until the border.

# ifdef REQUIRE
    Require("edge belongs to face", Has(elemX)) ;
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    Require("non-transposed mode", im == NORMAL_IM) ;
    InMethod("Quad::CopyAddValuesFromEdge(elemX,indx,indy,icx,icy,cam,im,nc") ;
# endif

    ElementaryField *eX, *eY, *eWork ;
    ParentElement   *borderParent, *mortarParent ;
    Matrix          *Q ;
    RealVector      *valX, *valY, *valWork, *valXX ;
    int             location, is, incr, np1, np2, dir, npX, npY, pX, pY ;
    int             foo1=cam, foo2=im ;// just to avoid a compilation warning

    eY  = GetElementaryField(indy) ;
    eX  = elemX->GetElementaryField(indx) ;

    pY  = procNumber ;
    pX  = elemX->GetProcessorNumber() ;

    npX = eX->GetParentElement(1)->GetNbCollocationPoints() ;
    np1 = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    np2 = eY->GetParentElement(2)->GetNbCollocationPoints() ;

    // 1. Determine how to address the 2D array of values

    location = SearchLocationOf(elemX) ;
    if (location == SOUTH) {
	is   = 1 ;
	incr = 1 ;
	dir  = 1 ;
	npY  = np1 ;
    }
    else if (location == NORTH) {
	is   = np1*(np2-1) + 1 ;
	incr = 1   ;
	dir  = 1 ;
	npY  = np1 ;
    }
    else if (location == WEST) {
	is   = 1 ;
	incr = np1 ;
	dir  = 2 ;
	npY  = np2 ;
    }
    else if (location == EAST) {
	is   = np1 ;
	incr = np1 ;
	dir  = 2 ;
	npY  = np2 ;
    }

    if (Mpi_rank == pY && rm == SIZE_RM) {
	Mpi_receiveBufferSizes[pX] += npY ;
	return ;
    }

    if (GetOrientationOfEdge(location) == -1)
	incr = -incr ;

    // 2. Do re-interpolation if required

    if (Mpi_rank == pY)                                  // for elemY's processor
	valY  = eY->GetComponent(icy) ;
    if (Mpi_rank == pX)                                  // for elemX's processor
	valX = eX->GetComponent(icx) ;

    if (eX->GetParentElement(1) == eY->GetParentElement(dir)) {
	// 2a. Conformal case: no re-interpolation

	valXX = valX ;
	eWork = NULL ;
    }

    else {
	// 2b. Nonconformal case: re-interpolate eX (1D) into eWork (1D); eWork has
	//     same interpolation as eY (2D) in the dir-th direction

	if (Mpi_rank == pX) {                              // for elemX's processor
	    eWork = eX->GetWork() ;
	    eWork->SetParentElement(eY->GetParentElement(dir),1) ;
	    valWork = eWork->GetComponent(1) ;
	    if (nc == MORTAR_NC) {
		// border = eWork, mortar = eX
		borderParent = eWork->GetParentElement(1) ;    // same parent as y
		mortarParent = eX->GetParentElement(1) ;       // same parent as x
		Q = elemX->GetMatrixQ(borderParent,mortarParent) ;
		Q->MatVec(valX,valWork) ;
	    }
	    else  // pointwise non-conformity
		eWork->CopyInterpolateFrom(eX,icx,1) ;
	    valXX = valWork ;
	}
    }

    // 3. Copy or add the values

    if (pX == pY) {                                   // same processor
	if (cam == COPY_CAM)
	    DCOPY(&npY, &valXX->At(1), &ione, &valY->At(is), &incr) ;
	else
	    DAXPY(&npY, &one, &valXX->At(1), &ione, &valY->At(is), &incr) ;
    }

    else {                                            // different processors
	if (Mpi_rank == pX)                             // for elemX's processor 
	    Mpi_SendToBuffer(npY, &valXX->At(1), ione, pY) ;
	else {                                          // for elemY's processor
	    if (cam == COPY_CAM)
		Mpi_ReceiveFromBuffer(npY, &valY->At(is), incr, pX, 0) ;
	    else
		Mpi_ReceiveFromBuffer(npY, &valY->At(is), incr, pX, 1) ;
	}
    }

    if (Mpi_rank == pX && eWork)
	ElementaryField::Retrieve(eWork) ;
}


void Quad :: CopyAddValuesFromFace (Face* elemX, 
                                    int indx, int indy, int icx, int icy,
                                    CopyAddMode cam, InterpMode im, 
                                    NonConform, ReceiveMode) 
{
    // Face from face version of method CopyAddValuesFrom.

# ifdef REQUIRE
    Require("same object", elemX == this) ;
    Require("correct processor", Mpi_rank == procNumber) ;
    InMethod("Quad::CopyAddValuesFromFace") ;
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


void Quad :: CopyAddValuesFromVolume (Volume* elemX, 
                                      int indx, int indy, int icx, int icy,
                                      CopyAddMode cam, InterpMode im, 
                                      NonConform nc, ReceiveMode rm) 
{
    // Face from volume version of method CopyAddValuesFrom.
    // No interpolation, extreme values are chosen. This means implicitly that
    // the collocation schemes of the 'indy'-th elementary field of the receiver 
    // are defined until the border.

# ifdef REQUIRE
    Require("face belongs to volume", elemX->Has(this)) ;
    Require("if mortar nonconformity then transposed interpolation", 
	    nc != MORTAR_NC || im == TRANSPOSED_IM) ;
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Quad:CopyAddValuesFromVolume(elemX,indx,indy,icx,icy,cam,in,nc)") ;
# endif

    ElementaryField *eY, *eX, *eWork1, *eWork2 ;
    RealVector      *valY, *valWork1, *valWork2 ;
    int             pY, dirKsi, dirEta, npKsi, npEta, southWest, southEast, 
	northWest, northEast ;
 
    pY = procNumber ;
    eY = GetElementaryField(indy) ;
    eX = elemX->GetElementaryField(indx) ;

    // 1. Determine which 2 directions out the 3 directions of y are relevant

    elemX->GetDirectionsOf(this,dirKsi,dirEta) ;

    // 2. Copy eX (3D) into a field eWork1 (2D) which has the same structure as
    //    eY

    if (eY->GetParentElement(1) == eX->GetParentElement(dirKsi) &&
	eY->GetParentElement(2) == eX->GetParentElement(dirEta)) {

	// 2.a: Conformal: x (3D) is contracted (including adequate rotations) into
	//      eWork1 (2D)

	eWork1 = elemX->Extract(indx,icx,this,rm) ;
	if (Mpi_rank == pY && rm == SIZE_RM)        // special case (handled by the 
	    return ;                                  // call to Extract hereabove)

	valWork1 = eWork1->GetComponent(1) ;
    }

    else {

	// 2.b: Non conformal:
	//   eX (3D) is contracted (including adequate rotations) into eWork2 (2D);
	//   eWork2 is projected into eWork1, with eWork1 having same structure as
	//   eY

	if (elemX->GetProcessorNumber() != procNumber)
	    Error("Face::CopyAddValuesFromVolume",
		  "non-conformity not implemented in parallel") ;

	eWork1   = eY->GetWork() ;
	valWork1 = eWork1->GetComponent(1) ;

	eWork2   = elemX->Extract(indx,icx,this,rm) ;
	valWork2 = eWork2->GetComponent(1) ;

	if (nc == MORTAR_NC)
	    ApplyQ(eWork2,eWork1,1,1,TRANSPOSED_IM) ;    // eWork1 = Q(transp).eWork2
	else {
	    if (im == TRANSPOSED_IM) 
		eWork1->CopyInterpolateTFrom(eWork2,1,1) ;
	    else
		eWork1->CopyInterpolateFrom(eWork2,1,1) ;
	}
	ElementaryField::Retrieve(eWork2) ;
    }

    // 3. Need to subtract 2/3 corner values since they are added from the three
    //    faces of the volume

    if (Mpi_rank == pY) {
	valY = eY->GetComponent(icy) ;
	if (cam == ADD_CAM) {
	    npKsi     = eY->GetParentElement(1)->GetNbCollocationPoints() ;
	    npEta     = eY->GetParentElement(2)->GetNbCollocationPoints() ;
	    southWest = 1 ;
	    southEast = npKsi ;
	    northWest = npKsi*(npEta-1) + 1 ;
	    northEast = npKsi*npEta ;
	    valY->At(southWest) -= TWO_THIRD * valWork1->At(southWest) ;
	    valY->At(southEast) -= TWO_THIRD * valWork1->At(southEast) ;
	    valY->At(northWest) -= TWO_THIRD * valWork1->At(northWest) ;
	    valY->At(northEast) -= TWO_THIRD * valWork1->At(northEast) ;
	}
    }

    // 4. Copy/add work1 (2D) into y (2D)

    if (Mpi_rank == pY) {
	if (cam == ADD_CAM)
	    valY->Add(valWork1) ;
	else
	    valY->CopyFrom(valWork1) ;
    }
    ElementaryField::Retrieve(eWork1) ;
}


Point* Quad :: Eval (real r, real s)
{
    // Returns a new point P(r,s), with 'r' and 's' both in [-1,1].

# ifdef REQUIRE
    Require("r and s are in [-1,1]", r>=-1. && r<=1. && s>=-1. && s<=1.) ;
    InMethod("Quad::Eval(r,s)") ;
# endif
  
    Point *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8 ;
    real   p1r, p2r, p1s, p2s, sWest, sEast, rSouth, rNorth ;

    // Interpolation functions at r and s
    p1r = Phi1(r) ;
    p2r = Phi2(r) ;
    p1s = Phi1(s) ;
    p2s = Phi2(s) ;

    // Contribution of the edges (use -s instead of s if negative orientation)
    sWest  = orientations->At(WEST) ==+1 ? s:-s ;
    sEast  = orientations->At(EAST) ==+1 ? s:-s ;
    rSouth = orientations->At(SOUTH)==+1 ? r:-r ;
    rNorth = orientations->At(NORTH)==+1 ? r:-r ;

    p1 = GetEdge(WEST) ->Eval(sWest) ->Multiply(p1r) ;
    p2 = GetEdge(EAST) ->Eval(sEast) ->Multiply(p2r) ;
    p3 = GetEdge(SOUTH)->Eval(rSouth)->Multiply(p1s) ;
    p4 = GetEdge(NORTH)->Eval(rNorth)->Multiply(p2s) ;

    // Contribution of the vortices
    p5 = GetVertex(1)->GetPoint()->Times(-p1r*p1s) ;       // South West
    p6 = GetVertex(4)->GetPoint()->Times(-p1r*p2s) ;       // North West
    p7 = GetVertex(2)->GetPoint()->Times(-p2r*p1s) ;       // South East
    p8 = GetVertex(3)->GetPoint()->Times(-p2r*p2s) ;       // North East

    // Sum the contributions
    p1->Add(p2)->Add(p3)->Add(p4)->Add(p5)->Add(p6)->Add(p7)->Add(p8) ;

    unref(p2) ;
    unref(p3) ;
    unref(p4) ;
    unref(p5) ;
    unref(p6) ;
    unref(p7) ;
    unref(p8) ;

    return p1 ;
}


Edge* Quad :: ExtractEdgeBetween (real r1, real r2, real s1, real s2)
{
    // Returns a new edge, the restriction of the receiver to [r1,r2]x[s1,s2],
    // where either r1=r2 or s1=s2.
    // If (r1,r2,s1,s2) actually defines
    //  - a portion of a contour edge of the receiver, returns the restriction of
    //    that contour (ie, a Line, or a Circle, etc); 
    //  - a portion of curve inside the receiver, returns a TransfiniteEdge.

# ifdef REQUIRE
    Require("r1,r2,s1,s2 range in [-1,1]",
	    r1>=-1. && r1<=1. && r2>=-1. && r2<=1. &&
	    s1>=-1. && s1<=1. && s2>=-1. && s2<=1.) ;
    Require("r1 <= r2 and s1 <= s2", r1 <= r2 && s1 <= s2) ;
    Require("r1 = r2 or s1 = s2", r1 == r2 || s1 == s2) ;
# endif

    Edge   *answer ;
    Vertex *v1, *v2 ;
    Point  *p1, *p2 ;

    if (r1 == r2) {                                          // "vertical"
	if (r1 == ONE_MINUS) {                                 // on the West edge
	    if (orientations->At(WEST) == +1)
		answer = GetEdge(WEST)->DuplicateBetween(s1,s2) ;
	    else
		answer = GetEdge(WEST)->DuplicateBetween(-s2,-s1) ;
	}
	else if (r1 == ONE) {                                  // on the East edge
	    if (orientations->At(EAST) == +1)
		answer = GetEdge(EAST)->DuplicateBetween(s1,s2) ;
	    else
		answer = GetEdge(EAST)->DuplicateBetween(-s2,-s1) ;
	}
	else {                                                 // in between
	    if (curvedInteriorEdges == true)
		answer = new TransfiniteEdge(this,1,r1,s1,s2) ;
	    else {
		p1 = Eval(r1,s1) ;
		p2 = Eval(r2,s2) ;
		v1 = new Vertex(p1) ;
		v2 = new Vertex(p2) ;
		answer = new Line(v1,v2) ;
		unref(p1) ;
		unref(p2) ;
		unref(v1) ;
		unref(v2) ;
	    }
	}
    }

    else {                                                   // "horizontal"
	if (s1 == ONE_MINUS) {                                 // on the South edge
	    if (orientations->At(SOUTH) == +1)
		answer = GetEdge(SOUTH)->DuplicateBetween(r1,r2) ;
	    else
		answer = GetEdge(SOUTH)->DuplicateBetween(-r2,-r1) ;
	}
	else if (s1 == ONE) {                                  // on the North edge
	    if (orientations->At(NORTH) == +1)
		answer = GetEdge(NORTH)->DuplicateBetween(r1,r2) ;
	    else
		answer = GetEdge(NORTH)->DuplicateBetween(-r2,-r1) ;
	}
	else {                                                 // in between
	    if (curvedInteriorEdges == true)
		answer = new TransfiniteEdge(this,2,s1,r1,r2) ;
	    else {
		p1 = Eval(r1,s1) ;
		p2 = Eval(r2,s2) ;
		v1 = new Vertex(p1) ;
		v2 = new Vertex(p2) ;
		answer = new Line(v1,v2) ;
		unref(p1) ;
		unref(p2) ;
		unref(v1) ;
		unref(v2) ;
	    }
	}
    }

    return answer ;
}


Mesh2D* Quad :: GenerateMesh () 
{
    // Using the two distributions of the receiver, returns a new mesh whose 
    // faces and vertices are generated by Gordon-Hall transfinite interpolations
    // between the receiver's edges.

    Vector<Edge*>   *eHoriz, *eVerti ;
    Vector<Vertex*> *v ;
    Mesh2D          *mesh ;
    Quad            *face ;
    Edge            *edge, *south, *north, *west, *east ;
    Distribution    *distrib1, *distrib2 ;
    Point           *newPoint ;
    real            r, s, r1, r2, s1, s2, absr ;
    int             i, j, idx, idxh, idxv, np1, np2, nseg1, nseg2 ;

    // Generation of the vertices

    distrib1 = GetDistribution(1) ;
    distrib2 = GetDistribution(2) ;
    np1      = distrib1->GetSize() ;
    np2      = distrib2->GetSize() ;
    v        = new Vector<Vertex*>(np1*np2) ;
    idx      = 1 ;
    for (j=1 ; j<=np2 ; j++) {
	s = distrib2->At(j) ;
	for (i=1 ; i<=np1 ; i++) {
	    r          = distrib1->At(i) ;
	    newPoint   = Eval(r,s) ;
	    v->At(idx) = new Vertex(newPoint) ;
	    unref(newPoint) ;
	    idx++ ;
	}
    }

    // Generation of the edges

    nseg1  = np1 - 1 ;
    nseg2  = np2 - 1 ;
    eHoriz = new Vector<Edge*>(np2*nseg1) ;
    eVerti = new Vector<Edge*>(np1*nseg2) ;
    idx    = 1 ;
    for (j=1 ; j<=np2 ; j++) {
	s    = distrib2->At(j) ;
	//    abss = fabs(s) ;
	for (i=1 ; i<=nseg1 ; i++) {
	    idxh = (j-1)*np1 + i ;
	    r1   = distrib1->At(i) ;
	    r2   = distrib1->At(i+1) ;
	    edge = ExtractEdgeBetween(r1,r2,s,s) ;
	    //   if (abss == ONE)                   // does not apply to transfinite edges
	    edge->SetVerticesTo(v->At(idxh),v->At(idxh+1)) ;
	    eHoriz->At(idx) = edge ;
	    idx++ ;
	}
    }
    idx = 1 ; 
    for (i=1 ; i<=np1 ; i++) {
	r    = distrib1->At(i) ;
	absr = fabs(r) ;
	for (j=1 ; j<=nseg2 ; j++) {
	    idxv = (j-1)*np1 + i ;
	    s1   = distrib2->At(j) ;
	    s2   = distrib2->At(j+1) ;
	    edge = ExtractEdgeBetween(r,r,s1,s2) ;
	    //  if (absr == ONE)                   // does not apply to transfinite edges
	    edge->SetVerticesTo(v->At(idxv),v->At(idxv+np1)) ;
	    eVerti->At(idx) = edge ;
	    idx++ ;
	}
    }

    // Generation of the faces
 
    mesh = new Mesh2D() ;
    for (j=1 ; j<=nseg2 ; j++)
	for (i=1 ; i<=nseg1 ; i++) {
	    idxh  = (j-1)*nseg1 + i ;
	    idxv  = (i-1)*nseg2 + j ;
	    south = eHoriz->At(idxh) ;
	    north = eHoriz->At(idxh+nseg1) ;
	    west  = eVerti->At(idxv) ;
	    east  = eVerti->At(idxv+nseg2) ;
	    face  = new Quad(south,north,west,east) ;
	    mesh->Put(face) ;
	    unref(face) ;
	}

    // cleaning
    for (i=1 ; i<=v->GetSize() ; i++)
	unref(v->At(i)) ;
    for (i=1 ; i<=eHoriz->GetSize() ; i++)
	unref(eHoriz->At(i)) ;
    for (i=1 ; i<=eVerti->GetSize() ; i++)
	unref(eVerti->At(i)) ;
    delete v ;
    delete eHoriz ;
    delete eVerti ;

    return mesh ;
}


int Quad :: GetDirectionOf (Edge* edge)
{
    // Returns the direction (1 or 2) of 'edge' in the receiver.

# ifdef REQUIRE
    Require("edge belongs to face", Has(edge)) ;
    InMethod("Quad::GetDirectionOf(edge)") ;
# endif

    int i ;

    i = edges->GetIndexOf(edge) ;
    if (i==SOUTH || i==NORTH)
	return 1 ;
    else
	return 2 ;
}


int Quad :: GetOrientationOfEdge (int i)
{
    // Returns +1 if the i-th edge is oriented positively, else returns -1.
    // See the sign convention for every edge in the header file.

# ifdef REQUIRE
    Require("'i' is in [1..4]", i >= 1 && i <= 4) ;
    InMethod("Quad::GetOrientationOfEdge(i)") ;
# endif

    return orientations->At(i) ;
}


void Quad :: GetSetInteriorDof (int index, int iop, int idof, real& val)
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
    InMethod("Hexa::GetSetInteriorDof(index,iop,idof,val)") ;
# endif

    ElementaryField *eY ;
    RealVector      *valY ;
    int             i, j, ic, np1, npi1, npi2, size, loc ;

    eY   = GetElementaryField(index) ;
    np1  = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    npi1 = eY->GetParentElement(1)->GetNbInteriorCollocationPoints() ;
    npi2 = eY->GetParentElement(2)->GetNbInteriorCollocationPoints() ;

    size = npi1 * npi2 ;
    ic   = (idof-1)/size + 1 ;
    idof = (idof-1)%size + 1 ;
    j    = (idof-1)/npi1 + 1 ;                        // j-th interior row
    i    = (idof-1)%npi1 + 1 ;                        // i-th interior column
    loc  = j*np1 + i + 1 ;

    valY = eY->GetComponent(ic) ;
    if (iop)
	valY->At(loc) = val ;
    else
	val = valY->At(loc) ;
}


Vertex* Quad :: GetVertex (int i)
{
    // Returns the i-th vertex of the receiver, numbered from 1 (South West) to
    // 4 (North West) counterclockwise.

# ifdef REQUIRE
    Require("'i' is in [1..4]", i >= 1 && i <= 4) ;
    InMethod("Quad::GetVertex(i)") ;
# endif

    int k, edge ;

    if (i == 1) {
	edge = SOUTH ;
	k    = +1 ;
    }
    else if (i == 2) {
	edge = SOUTH ;
	k    = -1 ;
    }
    else if (i == 3) {
	edge = NORTH ;
	k    = -1 ;
    }
    else {
	edge = NORTH ;
	k    = +1 ;
    }  

    if (k*orientations->At(edge) == +1)
	return edges->At(edge)->GetVertex1() ;
    else
	return edges->At(edge)->GetVertex2() ;
}


real Quad :: GetTypicalLength (int dir)
{
    // Returns a rough approximation of the mean length of the receiver in the
    // local 'dir'-th direction.
    // 
    //   THE ALGORITHM CAN BE SIMPLIFIED, SEE Hexa::GetTypicalLength(dir).

# ifdef REQUIRE
    Require("valid argument 'dir'", dir == 1 || dir == 2) ;
    InMethod("Quad::GetTypicalLength(dir)") ;
# endif

    Point *edge1, *edge2, *p1, *p2, *unit ;
    real  answer ;
 
    if (dir == 1) {
	edge1  = GetVertex(2)->GetPoint()->Minus(GetVertex(1)->GetPoint()) ;
	edge2  = GetVertex(3)->GetPoint()->Minus(GetVertex(4)->GetPoint()) ;
	p1     = Eval(-ONE,ZERO) ;
	p2     = Eval( ONE,ZERO) ;
    }
    else {
	edge1  = GetVertex(4)->GetPoint()->Minus(GetVertex(1)->GetPoint()) ;
	edge2  = GetVertex(3)->GetPoint()->Minus(GetVertex(2)->GetPoint()) ;
	p1     = Eval(ZERO,-ONE) ;
	p2     = Eval(ZERO, ONE) ;
    }
   
    unit   = p2->Minus(p1)->Normalized() ;
    answer = (fabs(edge1->Dot(unit)) + fabs(edge2->Dot(unit))) * HALF ;

    delete edge1 ;
    delete edge2 ;
    unref(p1) ;
    unref(p2) ;
    delete unit ;

    return answer ;
}


boolean Quad :: HasConsistentEdges ()
{
    // Returns True if the edges of the receiver are properly connected and if
    // they form a convex polygon.
    // Includes determining their orientation.

    Edge    *south, *north, *west, *east ;
    Vertex  *sw, *se, *nw, *ne ; 
    boolean    answer ;

    delete orientations ;                      // usually is already NULL
    orientations = new Vector<int>(4) ;
    orientations->SetValues(0) ;

    south = edges->At(SOUTH);
    north = edges->At(NORTH);
    west  = edges->At(WEST);
    east  = edges->At(EAST);

    // check connectivities and set orientations
    if (west->Has(south->GetVertex1()) && east->Has(south->GetVertex2())) {
	sw = south->GetVertex1() ;
	se = south->GetVertex2() ;
	orientations->At(SOUTH) = +1 ; }
    if (west->Has(south->GetVertex2()) && east->Has(south->GetVertex1())) {
	sw = south->GetVertex2() ;
	se = south->GetVertex1() ;
	orientations->At(SOUTH) = -1 ; }
    if (west->Has(north->GetVertex1()) && east->Has(north->GetVertex2())) {
	nw = north->GetVertex1() ;
	ne = north->GetVertex2() ;
	orientations->At(NORTH) = +1 ; }
    if (west->Has(north->GetVertex2()) && east->Has(north->GetVertex1())) {
	nw = north->GetVertex2() ;
	ne = north->GetVertex1() ;
	orientations->At(NORTH) = -1 ; }
    if (south->Has(west->GetVertex1()) && north->Has(west->GetVertex2()))
	orientations->At(WEST) = +1 ;
    if (south->Has(west->GetVertex2()) && north->Has(west->GetVertex1()))
	orientations->At(WEST) = -1 ;
    if (south->Has(east->GetVertex1()) && north->Has(east->GetVertex2()))
	orientations->At(EAST) = +1 ;
    if (south->Has(east->GetVertex2()) && north->Has(east->GetVertex1()))
	orientations->At(EAST) = -1 ;

    answer = orientations->At(SOUTH) && orientations->At(NORTH) && 
	orientations->At(WEST)  && orientations->At(EAST) ;

    if (answer == false ){
	for (int i=1; i<=4; i++)
	    printf("orientations de egde %d is %d \n",i,orientations->At(i));
	printf("SOUTH = %d WEST =%d \n",SOUTH,WEST );   
    }
    return answer ;
}


boolean Quad :: IsIdenticalTo (Element* elem)
{
    // Returns True if the receiver is geometrically identical to 'elem'.

    Face *f ;
    Quad *q ;
    Edge *south1, *north1, *west1, *east1, *south2, *north2, *west2, *east2 ;

    if (! elem->IsFace())
	return false ;
    else
	f = (Face*) elem ;

    if (f == this)
	return true ;
    if (! f->IsQuad())
	return false ;

    q      = (Quad*)f ;
    south1 = GetEdge(SOUTH) ;
    north1 = GetEdge(NORTH) ;
    west1  = GetEdge(WEST)  ;
    east1  = GetEdge(EAST)  ;
    south2 = q->GetEdge(SOUTH) ;
    north2 = q->GetEdge(NORTH) ;
    west2  = q->GetEdge(WEST)  ;
    east2  = q->GetEdge(EAST)  ;

    if (south2->IsIdenticalTo(south1)) 
	return ((west2->IsIdenticalTo(west1) && east2->IsIdenticalTo(east1)) ||
		(west2->IsIdenticalTo(east1) && east2->IsIdenticalTo(west1)))
	    && north2->IsIdenticalTo(north1) ;
    else if (south2->IsIdenticalTo(east1))
	return ((west2->IsIdenticalTo(south1) && east2->IsIdenticalTo(north1)) ||
		(west2->IsIdenticalTo(north1) && east2->IsIdenticalTo(south1)))
	    && north2->IsIdenticalTo(west1) ;
    else if (south2->IsIdenticalTo(north1))
	return ((west2->IsIdenticalTo(east1) && east2->IsIdenticalTo(west1)) ||
		(west2->IsIdenticalTo(west1) && east2->IsIdenticalTo(east1)))
	    && north2->IsIdenticalTo(south1) ;
    else if (south2->IsIdenticalTo(west1))
	return ((west2->IsIdenticalTo(north1) && east2->IsIdenticalTo(south1)) ||
		(west2->IsIdenticalTo(south1) && east2->IsIdenticalTo(north1)))
	    && north2->IsIdenticalTo(east1) ;
    else
	return false ;
}


boolean Quad :: IsQuad ()
{
    // Returns True because the receiver is a Quad.

    return true ;
}


void Quad :: Print ()
{
    // Prints the receiver on standard output.

    Point *bary ;
    int   i ;

    printf("Quad with orientations:  ") ;
    if (orientations)
	for (i=1 ; i<=orientations->GetSize() ; i++)
	    printf("%+d ",orientations->At(i)) ;
    else
	printf("NULL") ;

    bary = GetBarycenter() ;
    printf("Barycenter: "); bary->Print() ;
    unref(bary) ;

    printf ("\nEdges:\n") ;
    for (i=1 ; i<=edges->GetSize() ; i++) {
	printf("edge %d:\n ",i) ;
	edges->At(i)->Print() ;
    }
}


void Quad :: SetCurvedInteriorEdges (boolean choice)
{
    // A class method.
    // Sets to 'choice' the class attribute 'curvedInteriorEdges' (see 
    // description in header file).

    curvedInteriorEdges = choice ;
}


void Quad :: SetToCoordinates (int index) 
{
    // Initializes the elementary field 'index' to the coordinates of the
    // receiver.
    // The receiver initializes as many components as it finds.
    // In every component the coordinates are stored according to the convention
    // given in class ElementaryField, i.e., row by row: first the South-most row
    // (from West to East), next the next row, etc, until the North-most row.

# ifdef REQUIRE
    Require("the field has at most 3 components", 
	    GetElementaryField(index)->GetNbComponents() <= 3) ;
    InMethod("Quad::SetToCoordinates(index)") ;
# endif

    ElementaryField *eY ;
    Point           *point ;
    RealVector      *collocationPoints1, *collocationPoints2 ;
    real            r, s ;
    int             i, j, icomp, np1, np2, ncomp, pos ;

    if (Mpi_rank != procNumber)
	return ;

    eY                 = GetElementaryField(index) ;
    ncomp              = eY->GetNbComponents() ;
    collocationPoints1 = eY->GetParentElement(1)->GetCollocationPoints();
    collocationPoints2 = eY->GetParentElement(2)->GetCollocationPoints();
    np1                = collocationPoints1->GetSize() ;
    np2                = collocationPoints2->GetSize() ;
    pos                = 0 ;
    for (j=1 ; j<=np2 ; j++) {
	s = collocationPoints2->At(j) ;
	for (i=1 ; i<=np1 ; i++) {
	    r     = collocationPoints1->At(i) ;
	    point = Eval(r,s) ;
	    pos++ ;
	    for (icomp=1; icomp<=ncomp ; icomp++) 
		eY->GetComponent(icomp)->At(pos) = point->GetCoordinate(icomp) ;
	    unref(point) ;
	}
    }
}

