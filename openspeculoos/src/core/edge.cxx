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

// edge.cxx

#include "core/edge.hxx"
using namespace std ;

static real one  = ONE ;
static int  ione = IONE ;


//____________________________________ Edge ___________________________________


Edge :: Edge ()
    : Element()
{
    // Empty constructor.

    vertex1 = NULL ;
    vertex2 = NULL ;
    InitializeDistributions() ;

    operatorsQ = new Vector<OperatorQ1D*>() ;
}


Edge :: Edge (Vertex* v1, Vertex* v2)
    : Element()
{
    // Constructor. Initializes the receiver to an edge with 'v1' and 'v2' as end
    // vertices.

    ref(v1) ;
    ref(v2) ;
    vertex1 = v1 ;
    vertex2 = v2 ;
    InitializeDistributions() ;

    operatorsQ = new Vector<OperatorQ1D*>() ;
}


Edge :: ~Edge ()
{
    // Destructor.
    // Programming note: about the loop in this method, see Vector::~Vector().

    int i ;

    unref(vertex1) ;
    unref(vertex2) ;
    for (i=1 ; i<=operatorsQ->GetSize() ; i++)
	delete operatorsQ->At(i) ;
    delete operatorsQ ;
}


void Edge :: CopyAddValuesFromEdge (Edge* elemX,int indx, int indy, int icx, int icy,CopyAddMode cam, InterpMode im,NonConform, ReceiveMode) 
{
    // Edge-from-edge version of method CopyAddValuesFrom.
    // Implemented only in the normal interpolation mode and normal (i.e.,
    // non-mortar) nonconformity mode.

# ifdef REQUIRE
    Require("same object", elemX == this) ;
    Require("correct processor", Mpi_rank == procNumber) ;
    InMethod("Edge::CopyAddValuesFromEdge(elemX,indx,indy,icx,icy,cam,,,)") ;
# endif

    ElementaryField *eY, *eX, *eWork ;

    eY = GetElementaryField(indy) ;
    eX = elemX->GetElementaryField(indx) ;

    if (eY->HasSameInterpolationAs(eX)) {

	// 1. Conformal case: simple copy

	if (cam == ADD_CAM)
	    eY->GetComponent(icy)->Add(eX->GetComponent(icx)) ;
	else
	    eY->GetComponent(icy)->CopyFrom(eX->GetComponent(icx)) ;
    }

    else {

	// 2. Nonconformal case: re-interpolate

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


void Edge :: CopyAddValuesFromFace (Face* elemX,int indx, int indy, int icx, int icy,CopyAddMode cam, InterpMode im,NonConform nc, ReceiveMode rm) 
{
    // Edge-from-face version of method CopyAddValuesFrom.
    // No interpolation, extreme values are chosen. This means implicitly that
    // the collocation scheme of the 'indy'-th elementary field of the receiver 
    // is defined until the border.

# ifdef REQUIRE
    Require("'elemX' is a quad", elemX->IsQuad()) ;
    Require("edge belongs to face", elemX->Has(this)) ;
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Edge:CopyAddValuesFromFace(elemX,inx,indy,icx,icy,cam,in,nc)") ;
# endif

    ParentElement   *borderParent, *mortarParent ;
    ElementaryField *eY, *eX, *eWorkY, *eWorkX ;
    Matrix          *Q ;
    RealVector      *valY, *valX, *valWorkY, *valWorkX ;
    int             location, is, incr, np1, np2, dir, pX, pY, npX, npY ;
   
    eY  = GetElementaryField(indy) ;
    eX  = elemX->GetElementaryField(indx) ;
    npY = eY->GetParentElement(1)->GetNbCollocationPoints() ;
    np1 = eX->GetParentElement(1)->GetNbCollocationPoints() ;
    np2 = eX->GetParentElement(2)->GetNbCollocationPoints() ;
    pX  = elemX->GetProcessorNumber() ;
    pY  = procNumber ;

    if (Mpi_rank == pY && rm == SIZE_RM) {
	Mpi_receiveBufferSizes[pX] += npY ;
	return ;
    }

    // 1. Determine how to address the 2D array of values

    location = elemX->SearchLocationOf(this) ;
    if (location == SOUTH) {
	is   = 1 ;
	incr = 1 ; 
	dir  = 1 ;
	npX  = np1 ;
    }
    else if (location == NORTH) {
	is   = np1*(np2-1) + 1 ;
	incr = 1 ;
	dir  = 1 ;
	npX  = np1 ;
    }
    else if (location == WEST) {
	is   = 1 ;
	incr = np1 ;
	dir  = 2 ;
	npX  = np2 ;
    }
    else if (location == EAST) {
	is   = np1 ;
	incr = np1 ;
	dir  = 2 ;
	npX  = np2 ;
    }

    if (((Quad*)elemX)->GetOrientationOfEdge(location) == -1)
	incr = - incr ;

    // 2. Copy x (2D) into a work field workY (1D)

    if (Mpi_rank == pY) {                                // for elemY's processor
	eWorkY   = eY->GetWork() ;
	valWorkY = eWorkY->GetComponent(1) ;
	valY     = eY->GetComponent(icy) ;
    }
    if (Mpi_rank == pX)                                  // for elemX's processor 
	valX = eX->GetComponent(icx) ;

    // cout << eY->GetParentElement(1) << '\t' << eX->GetParentElement(dir) << endl ; 

    if (eY->GetParentElement(1) == eX->GetParentElement(dir)) {
	// 2a. Conformal case: simple copy
	//     From x (2D) to workY (1D)

	if (pX == pY)                                      // same processor
	    DCOPY(&npY, &valX->At(is), &incr, &valWorkY->At(1), &ione) ;
	else {                                             // different processors
	    if (Mpi_rank == pX)                              // for elemX's processor
		Mpi_SendToBuffer(npY, &valX->At(is), incr, pY) ; 
	    else                                             // for elemY's processor
		Mpi_ReceiveFromBuffer(npY, &valWorkY->At(1), ione, pX) ;
	}

	if (Mpi_rank == pY && cam == ADD_CAM) {
	    // need to subtract 1/2 corner values since later on they are added twice
	    // to the corner, i.e., from the two edges (of the receiver) which have
	    // that corner in common
	    valY->At(1)   -= HALF * valWorkY->At(1) ;
	    valY->At(npY) -= HALF * valWorkY->At(npY) ;
	}
    }

    else {
	// 2b. Nonconformal case: copy with pointwise or mortar reinterpolation,
	//     - from x (2D) conformly to workX (1D), then
	//     - from workX (1D) non-conformly to workY (1D).

	//     1. Do not copy values from a poorer element. This situation arises
	//        in the rare case of method Field::EnforceMortarConstraints().

	if (cam == COPY_CAM && npX < npY && nc == POINTWISE_NC
	    && FlatField::GetPolynomialChoice() == RICH_PC) {
	    if (Mpi_rank == pY)
		ElementaryField::Retrieve(eWorkY) ;
	    return ;
	}

	//     2. Copy x (2D) conformly to workX (1D)

	if (Mpi_rank == pY) {
	    eWorkX = eY->GetWork() ;
	    eWorkX->SetParentElement(eX->GetParentElement(dir),1) ; 
	    valWorkX = eWorkX->GetComponent(1) ;
	}

	if (pX == pY)                                      // same processor
	    DCOPY(&npX, &valX->At(is), &incr, &valWorkX->At(1), &ione) ;
	else {                                             // different processors
	    if (Mpi_rank == pX)                              // for elemX's processor
		Mpi_SendToBuffer(npX, &valX->At(is), incr, pY) ; 
	    else                                             // for elemY's processor
		Mpi_ReceiveFromBuffer(npX, &valWorkX->At(1), ione, pX) ;
	}

	if (Mpi_rank == pY && cam == ADD_CAM) {
	    // need to subtract 1/2 corner values since later on they are added twice
	    // to the corner, i.e., from the two edges (of the receiver) which have
	    // that corner in common
	    valY->At(1)   -= HALF * valWorkX->At(1) ;
	    valY->At(npY) -= HALF * valWorkX->At(npX) ;
	}

	//     3. Copy workX (1D) non-conformly to workY (1D)

	if (Mpi_rank == pY) {
	    if (nc == MORTAR_NC) {
		borderParent = eWorkX->GetParentElement(1) ;        // same parent as x
		mortarParent = eWorkY->GetParentElement(1) ;        // same parent as y
	
		// cout << " Computing the Q matrix " << endl ;
		Q = GetMatrixQ(borderParent,mortarParent) ;
		// cout << " Computing the Q matrix " << endl ;
		if (im == TRANSPOSED_IM)
		    Q->MatTransVec(valWorkX,valWorkY) ;
		else
		    Q->MatVec(valWorkX,valWorkY) ;
	    }
	    else {                     // pointwise non-conformity
		if (im == TRANSPOSED_IM) 
		    eWorkY->CopyInterpolateTFrom(eWorkX,1,1) ;
		else
		    eWorkY->CopyInterpolateFrom(eWorkX,1,1) ;
	    }
	    ElementaryField::Retrieve(eWorkX) ;
	}
    }

    // 3. Copy workY (1D) conformly to y (1D)

    if (Mpi_rank == pY) {
	if (cam == ADD_CAM)
	    valY->Add(valWorkY) ;
	else
	    valY->CopyFrom(valWorkY) ;
	ElementaryField::Retrieve(eWorkY) ;
    }
}


void Edge :: CopyAddValuesFromVertex (Vertex* elemX, 
                                      int indx, int indy, int icx, int icy,
                                      CopyAddMode cam, InterpMode, 
                                      NonConform, ReceiveMode rm) 
{
    // Edge-from-vertex version of method CopyAddValuesFrom.
    // No interpolation, extreme values are chosen. This means implicitly that
    // the collocation scheme of the 'indy'-th elementary field of the receiver 
    // is defined until the border.

# ifdef REQUIRE
    Require("vertex belongs to edge", Has(elemX)) ;
    Require("valid parent element type",
	    GetElementaryField(indy)->GetParentElement(1)->IsDefinedOnBorders());
    Require("correct processor", Mpi_IsOn(procNumber) ||
	    Mpi_IsOn(elemX->GetProcessorNumber())) ;
    InMethod("Edge::CopyAddValuesFromVertex") ;
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
    location = SearchLocationOf(elemX) ;
    if (location == 1)                             // left vertex
	ival = 1 ;
    else                                           // right vertex
	ival = eY->GetNbDofs() ;

    if (pX == pY) {                                // same processor
	if (cam == ADD_CAM) 
	    eY->GetComponent(icy)->At(ival) += eX->GetComponent(icx)->At(1) ;
	else
	    eY->GetComponent(icy)->At(ival)  = eX->GetComponent(icx)->At(1) ;
    }

    else {                                         // different processors
	if (pX == Mpi_rank)                          // run by elemX's processor
	    Mpi_SendToBuffer(1, &eX->GetComponent(icx)->At(1), 1, pY) ;
	else {                                       // run by elemY's processor
	    Mpi_ReceiveFromBuffer(1, &val[0], 1, pX) ;
	    if (cam == ADD_CAM)
		eY->GetComponent(icy)->At(ival) += val[0] ;
	    else
		eY->GetComponent(icy)->At(ival)  = val[0] ;
	}
    }
}


Mesh1D* Edge :: GenerateMesh () 
{
    // Generates the mesh (edges and vertices) resulting from the geometry and
    // the distribution of the receiver.

    Mesh1D       *mesh ;
    Edge         *edge ;
    Distribution *distrib ;
    int          i, nbseg ;

    // initializations
    distrib = GetDistribution(1) ;
    mesh    = new Mesh1D() ;

    // create the edges and vertices
    nbseg = distrib->GetSize() - 1;
    for (i=1 ; i<=nbseg ; i++) {
	edge = DuplicateBetween(distrib->At(i),distrib->At(i+1)) ;
	mesh->Put(edge) ;
	unref(edge) ;
    }

    return mesh ;
}


Point* Edge :: GetBarycenter ()
{
    // A debugging tool.
    // Returns a new point, the barycenter of the receiver.

    return this->Eval(ZERO) ;
}


FullMatrix* Edge :: GetMatrixQ (ParentElement* border, ParentElement* mortar)
{
    // Returns the matrix of the mortar operator Q to transfer values from a
    // 1D field u_mortar(r) to a 1D field u_border(r), both defined on the 
    // receiver but with different parent elements.
    //    u_border = matrix_Q * u_mortar ,
    // where u_border and u_mortar stand for the values of u_border(r) and
    // u_mortar(r).
    //
    // The mortar condition is given by the 3 conditions
    //   . integral { psi(r) * (u_border(r) - u_mortar(r)) } = 0,
    //   . u_border (corner1) = u_mortar(corner1),
    //   . u_border (corner2) = u_mortar(corner2),
    // which can rewritten as
    //   . integral { psi(r) * u_border(r) } = integral { psi(r) * u_mortar(r) }
    //   . + the above corner conditions
    // or, matricially, as
    //   H_border u_border = H_mortar u_mortar
    // and thus
    //   Q = inv(H_border) * H_mortar
    //   
    // Q is computed if it has not been constructed yet for this border/mortar
    // pair.
    //
    // The jacobian is assumed to be constant, i.e., the receiver is assumed to
    // be rectilinear. Therefore this method introduces an error if the receiver 
    // is curved. A possible solution would be to add an argument 
    // "index_of_jacobian" in Element::CopyAddValuesFrom(...).

# ifdef REQUIRE
    Require("valid argument 'border'", border->IsDefinedOnBorders()) ;
    Require("valid argument 'mortar'", mortar->IsDefinedOnBorders()) ;
    InMethod("Edge::GetMatrixQ(border,mortar)") ;
# endif

    OperatorQ1D   *operatorQ ;
    ParentElement *psi, *integ ;
    Matrix        *H_mortar_to_integ, *H_psi_to_integ, *H_psi_to_border ;
    FullMatrix    *H_border_no_cc, *H_mortar_no_cc, *H_border, *H_mortar, *Q,
	*temp ;
    RealVector    *weights ;
    int           i, j, n, np_border, np_border_interior, np_mortar, np_integ ;

    // Search the Q matrix for this border/mortar pair

    n = operatorsQ->GetSize() ;
    for (i=1 ; i<=n ; i++) {
	operatorQ = operatorsQ->At(i) ;
	if (operatorQ->GetBorderParent() == border &&
	    operatorQ->GetMortarParent() == mortar)
	    return operatorQ->GetMatrix() ;
    }

    cout << " Computing the Q matrix " << endl ;

    // Q operator not found in the list, so construct it

    np_border = border->GetNbCollocationPoints() ;
    np_mortar = mortar->GetNbCollocationPoints() ;

    Q = new FullMatrix(np_border,np_mortar) ;

    // psi is the set of polynomials that vanish at the interior GLL collocation
    // points of 'border', i.e., the polynomial degree of psi is set to that of
    // 'border' - 2

    np_border_interior = np_border - 2 ;
    if (np_border_interior > 0)
	psi = ParentElementIGLL::GetParentElementOfDegree(np_border_interior-1) ;
    else if (np_border_interior == 0)
	psi = NULL ;
    else
	Error("Edge::GetMatrixQ(border,mortar)","negative polynomial degree") ;

    // 1. Compute H_border_no_cc = H_border without the corner conditions: 
    //      H_border_no_cc * u_border = integral { psi(r) u_border(r) }.
    //    The integration rule is 'border'.
    //    H_border_no_cc = weights_border * H_psi_to_border(T).

    if (np_border_interior > 0) {
	H_border_no_cc  = new FullMatrix(np_border_interior,np_border) ;
	H_psi_to_border = psi->GetInterpMatrixTo(border) ;
	weights         = border->GetWeights() ;
	for (i=1 ; i <= np_border_interior ; i++) 
	    for (j=1 ; j <= np_border ; j++) 
		H_border_no_cc->At(i,j) = H_psi_to_border->At(j,i) * weights->At(j) ;
    }
    else
	H_border_no_cc = NULL ;

    // 2. Compute H_mortar_no_cc = H_mortar without the corner conditions:
    //      H_mortar_no_cc * u_mortar = integral { psi(r) u_mortar(r) }.
    //    The integration rule is whichever of 'border' and 'mortar' is richer.
    //    H_mortar_no_cc = weights_integ * H_psi_to_integ(T) * H_mortar_to_integ.

    if (np_mortar > np_border) {
	integ    = mortar ;
	np_integ = np_mortar ;
    }
    else {
	integ    = border ;
	np_integ = np_border ;
    }

    if (np_border_interior > 0) {
	// a) temp = weights * H_psi_to_integ(T)
	temp           = new FullMatrix(np_border_interior,np_integ) ;
	H_psi_to_integ = psi->GetInterpMatrixTo(integ) ;
	weights        = integ->GetWeights() ;
	for (i=1 ; i <= np_border_interior ; i++) 
	    for (j=1 ; j <= np_integ ; j++) 
		temp->At(i,j) = H_psi_to_integ->At(j,i) * weights->At(j) ;
	// b) H_mortar_no_cc = temp * H_mortar_to_integ 
	H_mortar_no_cc    = new FullMatrix(np_border_interior,np_mortar) ;
	H_mortar_to_integ = mortar->GetInterpMatrixTo(integ) ;
	H_mortar_no_cc->MatMat(temp,H_mortar_to_integ) ;
	delete temp ;
    }
    else
	H_mortar_no_cc = NULL ;

    // 3. Add the corner conditions to get H_border and H_mortar

    H_border = new FullMatrix(np_border,np_border) ;
    H_mortar = new FullMatrix(np_border,np_mortar) ;
    for (j=1 ; j<=np_border ; j++) {
	H_border->At(1,j)         = ZERO ;
	H_border->At(np_border,j) = ZERO ;
	for (i=2 ; i<np_border ; i++)
	    H_border->At(i,j) = H_border_no_cc->At(i-1,j) ;
    }
    H_border->At(1,1)                 = ONE ;     
    H_border->At(np_border,np_border) = ONE ;     

    for (j=1 ; j<=np_mortar ; j++) {
	H_mortar->At(1,j)         = ZERO ;
	H_mortar->At(np_border,j) = ZERO ;
	for (i=2 ; i<np_border ; i++)
	    H_mortar->At(i,j) = H_mortar_no_cc->At(i-1,j) ;
    }
    H_mortar->At(1,1)                 = ONE ;
    H_mortar->At(np_border,np_mortar) = ONE ;

    printf("H_border: ") ; H_border->Print() ;
    printf("H_mortar: ") ; H_mortar->Print() ;

    // 4. Get Q by solving the square system H_border Q = H_mortar 

    H_border->Solve(Q,H_mortar) ;
  
    // 5. Store the solution Q and clean up

    operatorsQ->Put(new OperatorQ1D(border,mortar,Q)) ;

    delete H_border_no_cc ;
    delete H_mortar_no_cc ;
    delete H_border ;
    delete H_mortar ;

    return Q ;
}


Vertex* Edge :: GetOtherVertex (Vertex* v)
{
    // Returns the vertex (starting or ending) opposite to 'v'. 'v' is assumed to
    // be one of the two end vertices.

# ifdef REQUIRE
    Require("vertex 'v' belongs to the edge", v == vertex1 || v == vertex2) ;
    InMethod("Edge::GetOtherVertex(v)") ;
# endif

    if (v == vertex1)
	return vertex2 ;
    else
	return vertex1 ;
}


void Edge :: GetSetInteriorDof (int index, int iop, int idof, real &val)
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
    Require("valid parent element",
	    GetElementaryField(index)->GetParentElement(1)->IsDefinedOnBorders()) ;
    InMethod("Edge::GetSetInteriorDof(index,iop,idof,val)") ;
# endif

    int ic, size, loc ;

    size = GetElementaryField(index)->GetParentElement(1)
	->GetNbInteriorCollocationPoints() ;
    ic   = (idof-1)/size + 1 ;
    loc  = (idof-1)%size + 2 ;

    if (iop)
	GetElementaryField(index)->GetComponent(ic)->At(loc) = val ;
    else
	val = GetElementaryField(index)->GetComponent(ic)->At(loc) ;
}


Element* Edge :: GetSubelement (int i)
{
    // Returns the 'i'-th vertex of the receiver.

    return GetVertex(i) ;
}


Vertex* Edge :: GetVertex (int i)
{
 
    // Returns the i-th vertex of the receiver.

# ifdef REQUIRE
    Require("valid argument 'i'", i==1 || i==2) ;
    InMethod("Vertex::GetVertex(i)") ;
# endif

    if (i == 1)
	return vertex1 ;
    else
	return vertex2 ;
}


Vertex* Edge :: GetVertexAt (Point* p)
{
    // Returns the receiver's vertex which coincides with 'p'.

# ifdef REQUIRE
    Require("'p' is a vertex of the edge", HasVertexAt(p)) ;
    InMethod("Edge::GetVertexAt(p)") ;
# endif

    if (vertex1->GetPoint()->IsIdenticalTo(p))
	return vertex1 ;
    else
	return vertex2 ;
}


boolean Edge :: HasVertexAt (Point* p)
{
    // Returns True if 'p' coincides with a vertex of the receiver.

    return (p->IsIdenticalTo(vertex1->GetPoint()) ||
	    p->IsIdenticalTo(vertex2->GetPoint())) ;
}
  

boolean Edge :: HasSameOrientationAs (Edge* edge)
{
    // Returns True if the end vertices of the receiver and 'edge' are equal,
    // else returns False.

    return (vertex1->IsIdenticalTo(edge->vertex1) && 
	    vertex2->IsIdenticalTo(edge->vertex2)) ;
}


boolean Edge :: IsCircle ()
{
    // Returns false, because the receiver is not a circle.

    return false ;
}


boolean Edge :: IsLine ()
{
    // Returns False, because the receiver is not a line.

    return false ;
}


boolean Edge :: IsTransfiniteEdge ()
{
    // Returns False, because the receiver is not a transfinite edge.

    return false ;
}


void Edge :: Print ()
{
    // Prints the receiver on the standard output.

    printf("Edge\n") ;
    printf("    starts at: ") ; vertex1->Print() ;
    printf("    stops  at: ") ; vertex2->Print() ;
}


void Edge :: PrintPointers ()
{
    // A debugging tool.
    // Prints the receiver on the standard output (pointers only).

    printf("    Edge %x with vertices %x and %x\n",this,vertex1,vertex2) ;
}


int Edge :: SearchLocationOf (Element* v)
{
    // Returns the position number of the vertex 'v' in the receiver.
    // Returns 0 if 'v' does not belong to the receiver.

# ifdef REQUIRE
    Require("'v' is a vertex", v->IsVertex()) ;
    InMethod("Edge::SearchLocationOf(v)") ;
# endif

    if (v == vertex1)
	return 1 ;
    else if (v == vertex2)
	return 2 ;
    else
	return 0 ;
}


void Edge :: SetToCoordinates (int index) 
{
    // Sets the elementary field 'index' to the coordinates of the receiver.
    // The receiver initializes as many components as it finds.

# ifdef REQUIRE
    Require("field with valid number of components", 
	    GetElementaryField(index)->GetNbComponents() <= 3) ;
    InMethod("Edge::SetToCoordinates(index)") ;
# endif 

    ElementaryField *ef ;
    Point           *point ;
    RealVector      *collocationPoints ;
    real            r ;
    int             i, icomp, np, ncomp ;

    if (Mpi_rank != procNumber)
	return ;

    ef                = GetElementaryField(index) ;
    ncomp             = ef->GetNbComponents() ;
    collocationPoints = ef->GetParentElement(1)->GetCollocationPoints();
    np                = collocationPoints->GetSize() ;
    for (i=1 ; i<=np  ; i++) {
	r     = collocationPoints->At(i) ;
	point = Eval(r) ;
	for (icomp=1; icomp<=ncomp ; icomp++) 
	    ef->GetComponent(icomp)->At(i) = point->GetCoordinate(icomp) ;
	unref(point) ;
    }
}


void Edge :: SetVertex (int i, Vertex* v)
{
    // Sets the i-th vertex of the receiver to 'v'. 'v' must be geometrically
    // coincident with the vertex it replaces.

# ifdef REQUIRE
    Require("valid argument 'i'", i==1 || i==2) ;
    Require("'v' coincides with vertex", (i==1 && v->IsIdenticalTo(vertex1)) ||
	    (i==2 && v->IsIdenticalTo(vertex2))) ;
    InMethod("Edge::SetVertex(i,v)") ;
# endif

    ref(v) ;

    if (i == 1) {
	unref(vertex1) ;
	vertex1 = v ;
    }
    else {
	unref(vertex2) ;
	vertex2 = v ;
    }
}


void Edge :: SetVerticesTo (Vertex* v1, Vertex* v2)
{
    // Replaces vertex 1 with that of 'v1' and 'v2' that is geometrically coinci-
    // dent with vertex 1, and replaces vertex 2 with the other one.

# ifdef REQUIRE
    Require("v1 and v2 coincide with the vertices",
	    (v1->IsIdenticalTo(vertex1) && v2->IsIdenticalTo(vertex2)) ||
	    (v1->IsIdenticalTo(vertex2) && v2->IsIdenticalTo(vertex1))) ;
    InMethod("Edge::SetVerticesTo(v1,v2)") ;
# endif

    ref(v1) ;
    ref(v2) ;

    if (v1->IsIdenticalTo(vertex1)) {
	unref(vertex1) ;
	unref(vertex2) ;
	vertex1 = v1 ;
	vertex2 = v2 ;
    }
    else {
	unref(vertex1) ;
	unref(vertex2) ;
	vertex1 = v2 ;
	vertex2 = v1 ;
    }
}


void Edge :: Substitute (Vertex* source, Vertex* target)
{
    // In the definition of the receiver, substitutes 'source' by 'target',
    // where 'source' is one of the end vertices of the receiver.

# ifdef REQUIRE
    Require("has 'source'", source == vertex1 || source == vertex2) ;
    Require("'target' and 'source' coincide", target->IsIdenticalTo(vertex1) ||
	    target->IsIdenticalTo(vertex2)) ;
    InMethod("Edge::Substitute(source,target)") ;
# endif

    ref(target) ;
  
    if (source == vertex1) {
	unref(vertex1) ;
	vertex1 = target ;
    }
    else {
	unref(vertex2) ;
	vertex2 = target ;
    }
}


void Edge :: SubstituteCommonVertices (Edge* edge)
{
    // For each vertex v1 of the receiver that is identical to a vertex v2 of
    // 'edge', makes the receiver point to v2 rather than to v1.

    Vertex *v1, *v2 ;
    int    i, j ;

    for (i=1 ; i<=2 ; i++) {
	v2 = edge->GetVertex(i) ;
	for (j=1 ; j<=2 ; j++) {
	    v1 = GetVertex(j) ;
	    if (v1->IsIdenticalTo(v2))
		Substitute(v1,v2) ;
	}
    }
}


//___________________________________ Line ____________________________________


Line :: Line (Vertex* v1, Vertex* v2)
    : Edge (v1, v2)
{
    // Constructor. Initializes the receiver to a segment of line between 'v1'
    // and 'v2'.
}


Line :: ~Line ()
{
    // Destructor.
}


Edge* Line :: DuplicateBetween (real r1, real r2)
{
    // Returns a new straight edge, the restriction of the receiver to [r1,r2].
    // Remark: C++ is a shit; it should allow this method to return explicitly a
    //         Line, not just an Edge. Identity of signatures in base class Edge
    //         and in class Line prevents it.

# ifdef REQUIRE
    Require("'r1' < 'r2'", r1 < r2) ;
    Require("[r1,r2] is in [-1,1]", r1 >= ONE_MINUS || r2 <= ONE) ;
    InMethod("Line::DuplicateBetween(r1,r2)") ;
# endif

    Line   *answer ;
    Vertex *v1, *v2 ;
    Point  *p1, *p2 ;

    p1 = Eval(r1) ;
    p2 = Eval(r2) ;
    v1 = new Vertex(p1) ;
    v2 = new Vertex(p2) ;
    answer = new Line(v1,v2) ;

    unref(p1) ;
    unref(p2) ;
    unref(v1) ;
    unref(v2) ;

    return answer ;
}


Point* Line :: Eval (real r)
{
    // Returns a new point, the receiver evaluated at abscissa 'r'.
    // 'r' is usually in [-1,1].

    Point *p1, *p2 ;

    p1 = vertex1->GetPoint()->Times(Phi1(r)) ;
    p2 = vertex2->GetPoint()->Times(Phi2(r)) ;
    p1->Add(p2) ;

    unref(p2) ;
    return p1 ;
}


real Line :: GetLength ()
{
    // Returns the length of the receiver.

    Point *p1, *p2 ;
    real  proj, answer ;
    int   i ;

    p1     = vertex1->GetPoint() ;
    p2     = vertex2->GetPoint() ;
    answer = ZERO ;
    for (i=1 ; i<=3 ; i++) {
	proj    = p2->GetCoordinate(i) - p1->GetCoordinate(i) ;
	answer += proj * proj ;
    }

    return sqrt(answer) ;
}


real Line :: GetTypicalLength (int dir)
{
    // Returns the length of the receiver.

# ifdef REQUIRE
    Require("'dir' = 1", dir == 1) ;
    InMethod("Line::GetMeanDirection(dir)") ;
# endif

    int foo = dir ;                 // just to avoid a compilation warning!

    return GetLength() ;
}


boolean Line :: IsIdenticalTo (Element* elem)
{
    // Returns True if the receiver and 'elem' are identical (not necessarily the
    // same object).

    Edge *edge ;

    if (! elem->IsEdge())
	return false ;
    else
	edge = (Edge*) elem ;

    if (edge == this)
	return true ;
    else if (edge->IsLine())
	return (edge->HasVertexAt(vertex1->GetPoint()) &&
		edge->HasVertexAt(vertex2->GetPoint())) ;
    else
	return false ;
}


boolean Line :: IsLine ()
{
    // Returns True, because the receiver is a line.

    return true ;
}


void Line :: Print ()
{
    // Prints the receiver on standard output.

    printf ("  Line: ") ;
    Edge::Print() ;
}


//__________________________________ Circle ___________________________________

Circle :: Circle (Vertex* v1, Vertex* v2, Point* c)
    : Edge (v1, v2)
{
    // Constructor. Initializes the receiver to a circular arc starting at 'v1'
    // and ending at 'v2', with center a copy of 'c'.
    // 'v1', 'v2' and 'c' must lie in the same horizontal plane Z=cst.

    Point *tmp1, *tmp2 ;
    real  r1 ;

    center = c->Duplicate() ;
    tmp1   = v1->GetPoint()->Minus(c) ;
    tmp2   = v2->GetPoint()->Minus(c) ;
    r1     = tmp1->GetNorm() ;
    radius = r1 ;
    angle1 = atan2(tmp1->GetY(), tmp1->GetX()) ;
    if (tmp1->GetY()<0.0) angle1 = angle1 + 2.0*PI ;
    angle2 = atan2(tmp2->GetY(), tmp2->GetX()) ;
    if (tmp2->GetY()<0.0) angle2 = angle2 + 2.0*PI ;
  
    if (fabs(angle2-angle1)>PI){
	printf("\n WARNING: arc angle greater than PI! \n");
	if (angle1<1.0e-10) angle1 = angle1 + 2.0*PI;
	if (angle2<1.0e-10) angle2 = angle2 + 2.0*PI;


    }

# ifdef CHECK
    Point *tmp3 = v2->GetPoint()->Minus(v1->GetPoint()) ;
    real r2     = tmp2->GetNorm() ;
    real d      = tmp3->GetNorm() ;
    if (fabs(r1-r2)/r2 > DEFAULT_TOLERANCE || d < r1/100)
	Error("Circle::Circle","the points are not defining a circle") ;
    if (fabs(v1->GetPoint()->GetZ()-v2->GetPoint()->GetZ()) > DEFAULT_TOLERANCE ||
	fabs(v1->GetPoint()->GetZ()-c->GetZ()) > DEFAULT_TOLERANCE)
	Error("Circle::Circle","v1,v2,c not in the same Z-plane") ;
    unref(tmp3) ;
# endif

    unref(tmp1) ;
    unref(tmp2) ;

}

Circle :: ~Circle ()
{
    // Destructor.

    unref(center) ;
}

Edge* Circle :: DuplicateBetween (real r1, real r2)
{
    // Returns a new circular arc, the restriction of the receiver to [r1,r2],
    // with 'r1' and 'r2' both in [-1,1].

# ifdef REQUIRE
    Require("r1 < r2", r1 < r2) ;
    InMethod("Circle::DuplicateBetween(r1,r2)") ;
# endif

    Circle *answer ;
    Vertex *v1, *v2 ;
    Point  *p1, *p2 ;

# ifdef CHECK
    if (r1 < -1. || r2 > 1.)
	Warning("Line::DuplicateBetween(r1,r2)",
		"r1 and r2 expected to be in [-1,1]") ;
# endif

    p1 = Eval(r1) ;
    p2 = Eval(r2) ;
    v1 = new Vertex(p1) ;
    v2 = new Vertex(p2) ;
    answer = new Circle(v1,v2,center) ;

    unref(p1) ;
    unref(p2) ;
    unref(v1) ;
    unref(v2) ;

    return answer ;
}


Point* Circle :: Eval (real r)
{
    // Returns a new point, the receiver evaluated at abscissa 'r', with 'r' 
    // typically in [-1,1].
    // r=-1 -> start point.
    // r= 1 -> stop point.

    Point *answer ;
    real  theta ;

    theta  = Phi1(r)*angle1 + Phi2(r)*angle2 ;
    answer = new Point(radius*cos(theta),radius*sin(theta),ZERO) ;
    answer->Add(center) ;

    return answer ;
}


boolean Circle :: IsCircle ()
{
    // Returns True, because the receiver is a circle.

    return true ;
}


boolean Circle :: IsIdenticalTo (Element* elem)
{
    // Returns True if the receiver and 'elem' are identical (not necessarily the
    // same object).

    Edge   *edge ;
    Circle *circle ;

    if (! elem->IsEdge())
	return false ;
    else
	edge = (Edge*)elem ;

    if (edge== this)
	return true ;
    else if (edge->IsCircle()) {
	circle = (Circle*)edge ;
	return circle->HasVertexAt(vertex1->GetPoint()) &&
	    circle->HasVertexAt(vertex2->GetPoint()) &&
	    circle->GetCenter()->IsIdenticalTo(center) ;
    }
    else
	return false ;
}


void Circle :: Print ()
{
    // Prints the receiver on standard output.

    printf("  Circle: ") ; 
    Edge::Print() ;
    printf("    radius %6.2f   angle1 %6.1f degrees   angle2 %6.1f degrees\n",
	   radius,angle1*180/PI,angle2*180/PI);
    printf("    center: ") ; center->Print() ;
}


//___________________________ TransfiniteEdge _________________________________


TransfiniteEdge :: TransfiniteEdge (Quad* q, int dir, real RS0, 
                                    real A1, real A2)
    : Edge ()
{
    // Constructor. Initializes the receiver to a transfinite edge defined as the
    // restriction of 'q' (which defines a function F(r,s)) to:
    //   . F(r0,s), if 'dir'=1; r0 is given by 'RS0'; s ranges in [A1,A2], or
    //   . F(r,s0), if 'dir'=2. s0 is given by 'RS0'; r ranges in [A1,A2].

# ifdef REQUIRE
    Require("'dir' is 1 or 2",  dir == 1 || dir == 2) ;
    Require("RS0 is in [-1,1]", RS0 >= -1. && RS0 <= 1.) ;
    Require("A1 is in [-1,1]",  A1 >= -1. && A1 <= 1.) ;
    Require("A2 is in [-1,1]",  A2 >= -1. && A2 <= 1.) ;
    InMethod("TransfiniteEdge::TransfiniteEdge(q,dir,RS0,A1,A2)") ;
# endif

    Point *p1, *p2 ;

    quad             = q ;
    blockedDirection = dir ;
    rs0              = RS0 ;
    a1               = A1 ;
    a2               = A2 ;

    if (dir == 1) {
	p1 = quad->Eval(rs0,a1) ;
	p2 = quad->Eval(rs0,a2) ;
    }
    else {
	p1 = quad->Eval(a1,rs0) ;
	p2 = quad->Eval(a2,rs0) ;
    }

    vertex1 = new Vertex(p1) ;
    vertex2 = new Vertex(p2) ;
    unref(p1) ;
    unref(p2) ;
}


Edge* TransfiniteEdge :: DuplicateBetween (real , real )
{
    // Issues an error message.

    NotImplemented("TransfiniteEdge::DuplicateBetween(s1,s2)") ;
    return NULL ;
}


Point* TransfiniteEdge :: Eval (real rs)
{
    // Returns a new point, the receiver evaluated at 'rs', where 'rs' is the
    // free abscissa of the receiver (ie, s if blockedDirection=1, r else).
    // 'rs' should be in [-1,1].

    if (blockedDirection == 1)
	return quad->Eval(rs0,ParentToPhysic(a1,a2,rs)) ;
    else
	return quad->Eval(ParentToPhysic(a1,a2,rs),rs0) ;
}


boolean TransfiniteEdge :: IsIdenticalTo (Element* elem)
{
    // Returns True if the receiver and 'elem' are identical (not necessarily the
    // same object).
    // Exceptionally, pointer equality is checked here.

    return elem == this ;
}


boolean TransfiniteEdge :: IsTransfiniteEdge ()
{
    // Returns True, because the receiver is a transfinite edge.

    return true ;
}


void TransfiniteEdge :: Print ()
{
    // Prints the receiver on standard output.

    printf("  TransfiniteEdge on Quad %d: ",quad) ;
    Edge::Print() ;
    if (blockedDirection == 1)
	printf("    r0 = %6.3f   s1 = %6.3f   s2 = %6.3f\n",rs0,a1,a2) ;
    else
	printf("    s0 = %6.3f   r1 = %6.3f   r2 = %6.3f\n",rs0,a1,a2) ;
}


