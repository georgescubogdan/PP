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

// parent.cxx

#include "core/parent.hxx"
#include "core/listmat.hxx"
#include "core/matrix.hxx"
#include "core/jacobi.hxx"
#include "core/options.hxx"
#include <stdio.h>

#include <math.h>             // for the debugging methods only
#include <stdlib.h>           //             id.
#include <sstream>

const int MAX_DEGREE = 100 ;


//_______________________________ ParentElement _______________________________


ParentElement :: ParentElement (int n)
{
    // Constructor. Initializes the receiver to a parent element of degree 'n'.

    polynomialDegree        = n ;
    collocationPoints       = NULL ;
    quadratureWeights       = NULL ;
    derivMatrix             = NULL ;
    listInterpMatrices      = new ListMatrices() ;
    listInterpDerivMatrices = new ListMatrices() ;
}


ParentElement :: ~ParentElement ()
{
    // Destructor.

    delete collocationPoints ;
    delete quadratureWeights ;
    delete derivMatrix ;
    delete listInterpMatrices ;
    delete listInterpDerivMatrices ;
}


FullMatrix* ParentElement :: ComputeInterpDerivMatrixTo (ParentElement* p)
{
    // A debugging tool.
    // Computes (without resorting to the MIT library) the derivative matrix of
    // the receiver, reinterpolated into 'p'.
    // Its size is: nb_points(p) x nb_points(receiver).
    // Assumes that the collocation points of both rules are correct.
    // If 'p' is the receiver, then this method should return the same matrix as
    // the one given by the MIT library, i.e., given by method GetDerivMatrix.

    FullMatrix *answer ;
    real       xi, val ;
    int        i, j, dim, dimP ;

    dim  = collocationPoints->GetSize() ;
    dimP = p->collocationPoints->GetSize() ;

    answer = new FullMatrix(dimP,dim) ;

    for (j=1 ; j<=dim ; j++)                // for all hj of the receiver
	for (i=1 ; i<=dimP ; i++) {           // for all points xi of 'p'
	    xi  = p->collocationPoints->At(i) ;
	    val = EvaluatePolynomialDerivative(j,xi) ;    // compute hj(xi) 
	    answer->At(i,j) = val ;
	}

    return answer ;
}


FullMatrix* ParentElement :: ComputeInterpMatrixTo (ParentElement* p)
{
    // A debugging tool.
    // Computes (without resorting to the MIT library) the interpolation matrix
    // from the receiver to 'p'.
    // Its size is: nb_points(p) x nb_points(receiver).
    // Assumes that the collocation points of both rules are correct.

    FullMatrix *answer ;
    real       xi, val ;
    int        i, j, dim, dimP ;

    dim  = collocationPoints->GetSize() ;
    dimP = p->collocationPoints->GetSize() ;

    answer = new FullMatrix(dimP,dim) ;

    for (j=1 ; j<=dim ; j++)                // for all hj of the receiver
	for (i=1 ; i<=dimP ; i++) {           // for all points xi of 'p'
	    xi  = p->collocationPoints->At(i) ;
	    val = EvaluatePolynomial(j,xi) ;    // compute hj(xi) 
	    answer->At(i,j) = val ;
	}

    return answer ;
}


real ParentElement :: EvaluatePolynomial (int j, real x)
{
    // A debugging tool.
    // Evaluates the 'j'-th polynomial of the receiver at 'x', x in [-1,1].

# ifdef REQUIRE
    Require("valid argument 'j'", j >= 1 && j <= collocationPoints->GetSize()) ;
    InMethod("ParentElement::EvaluatePolynomial(j,x)") ;
# endif

    real xj, numer, denom ;
    int  k ;

    xj = collocationPoints->At(j) ;

    numer = ONE ;
    denom = ONE ;
    for (k=1 ; k<=collocationPoints->GetSize() ; k++)
	if (k != j) {
	    numer *= x - collocationPoints->At(k) ;
	    denom *= xj - collocationPoints->At(k) ;
	}

    return numer/denom ;
} 
  

real ParentElement :: EvaluatePolynomialDerivative (int j, real s)
{
    // A debugging tool.
    // Evaluates the 'j'-th polynomial of the receiver at 'x', x in [-1,1].
    //
    // The algorithm proceeds in the following way:
    // Consider the example where the receiver is a 4-th order rule with 
    // collocation points (a,b,d,d,e), all in [-1,1]. Argument j is in {1,2,3,4}.
    // If j=2, the polynomial is h2(s) = (s-a)(s-c)(s-d)(s-e).
    //
    // 1) expands the polynomial.
    //    coalesces the first 2 factors into 1:
    //      h2(s) = (A0 + A1.s + A2.s^2) (s-d) (s-e)
    //    do it again:
    //      h2(s) = (B0 + B1.s + B2.s^2 + B3.s^3) (s-e)
    //    and again:
    //      h2(s) = (C0 + C1.s + C2.s^2 + C3.s^3 + C4.s^4)
    // 2) computes the derivative at s='s' (the argument of this method).
    //
    // Note that actually h2(s) must be divided by the denominator
    // (b-a)(b-c)(b-d)(b-e).
    //
    // A polynomial of degree, say 2, h(s) = (A0 + A1.s + A2.s^2) / denom 
    // is stored in a vector (A0,A1,A2,denom).

# ifdef REQUIRE
    Require("valid argument 'j'", j >= 1 && j <= collocationPoints->GetSize()) ;
    InMethod("ParentElement::EvaluatePolynomialDerivative(j,x)") ;
# endif

    Vector<RealVector*> *polynomial ;
    RealVector          *factor, *first, *second, *newFirst, *expanded ;
    real                s_i, s_j, c_i, answer, contrib, denom ;
    int                 i, ii, k, nbCoeff ;

    // 1. Create the polynomial in its factored form (a set of monomials)
  
    polynomial = new Vector<RealVector*>(polynomialDegree) ;
    s_j        = collocationPoints->At(j) ;
    for (i=ii=1 ; i<=polynomialDegree+1 ; i++)
	if (i != j) {                              // s_j (i.e., b) must be skipped
	    // factor = (s - s_i) / (s_j - s_i)
	    s_i    = collocationPoints->At(i) ;
	    factor = new RealVector(3) ;
	    factor->At(1) = -s_i ;
	    factor->At(2) = ONE ;
	    factor->At(3) = s_j - s_i ;
	    polynomial->At(ii) = factor ;
	    ii++ ;
	}

    // 2. Expand the polynomial, by multiplying recursively the first factor by
    //    the second one (the second is always of degree 1)

    for (i=1 ; i<polynomialDegree ; i++) {
	first    = polynomial->At(1) ;
	second   = polynomial->At(2) ;
	newFirst = new RealVector(first->GetSize()+1) ;
	nbCoeff  = i+1 ;         // nb coeffs A0, A1, A2, etc, in first factor

	newFirst->At(1) = first->At(1) * second->At(1) ;    // B0
	for (k=2 ; k<=nbCoeff ; k++)                        // B1, B2, ..., Bi
	    newFirst->At(k) = first->At(k-1) + first->At(k)*second->At(1) ;  
	newFirst->At(nbCoeff+1) = first->At(nbCoeff) ;      // Bi+1

	newFirst->At(nbCoeff+2) = first->At(nbCoeff+1)*second->At(3) ;
	// denom(1) * denom(2)

	polynomial->At(1) = newFirst ;          // replace first with newFirst
	polynomial->RemoveIndex(2) ;
    }

    // 3. Compute derivative of the expanded polynomial (C0 + C1.s + C2.s^2 +
    //    C3.s^3 + C4.s^4) / denom, at s = 's' (the argument of this method).

    expanded = polynomial->At(1) ;
    denom    = expanded->At(polynomialDegree+2) ;
    answer   = ZERO ;
    for (i=1 ; i<=polynomialDegree ; i++) {
	c_i = expanded->At(i+1) ;
	if (i == 1)                         // beware of s^0 if s=0. !
	    contrib = c_i / denom ;
	else
	    contrib = c_i * i * pow(s,i-1) / denom ;
	answer += contrib ;
    }

    return answer ;
}
    
  
void ParentElement :: DeleteAllInstances ()
{
    // Deletes all instances (ie, objects) of the derived classes. This method is
    // typically invoked at the end of the program's execution, for releasing all
    // memory.

    ParentElementGL  ::DeleteAllInstances() ;
    ParentElementGLL ::DeleteAllInstances() ;
    ParentElementGLLI ::DeleteAllInstances() ;
    ParentElementIGLL::DeleteAllInstances() ;
}


FullMatrix* ParentElement :: GetInterpDerivMatrixTo (ParentElement* p)
{
    // Returns the interpolated-derivative matrix D from the receiver to 'p'.
    // Let u denote the discrete values on the receiver of a function u(x); then
    // D can be used to calculate v = D u, where v is the discrete projection on
    // 'p' of the derivative of u(x), u(x) given by u.

    FullMatrix *d ;

    d = (FullMatrix*) listInterpDerivMatrices->Search(p) ;
    if (d == NULL) {
	d = GetInterpMatrixTo(p)->Times(GetDerivMatrix()) ;
	listInterpDerivMatrices->Store(p,d) ;
    }

    return d ;
}


ParentElement* ParentElement :: GetParentElementOfDegree_ (int n)
{
    // Looks into the list of already created parent for the parent of degree 
    // 'n' of the same family as the receiver. If such parent does not exist, 
    // creates it. 
    // For example, if the receiver is a GLL parent element and n=4, then this 
    // method returns the GLL4 parent element.
    // The client should not delete the returned parent.
    // Programming note: an underscore (_) has been added to the name of this
    // method, because a static method cannot be virtual.

    if (GetType() == GLL_PT)
	return ParentElementGLL::GetParentElementOfDegree(n) ;
    
    else if (GetType() == GLLI_PT)
	return ParentElementGLLI::GetParentElementOfDegree(n) ;  

    else if (GetType() == IGLL_PT)
	return ParentElementIGLL::GetParentElementOfDegree(n) ;

    else if (GetType() == GL_PT)
	return ParentElementGL::GetParentElementOfDegree(n) ;

    else {
	Error("ParentElement::GetParentElementOfDegree_(n)","unknown parent type");
	return NULL ;
    }
}


RealVector* ParentElement :: GetWeights ()
{
    // Returns the array containing the integration weight associated with every
    // collocation point of the receiver.

    return quadratureWeights ;
}


boolean ParentElement :: IsDefinedOnBorders ()
{
    // Returns False (default), since the receiver has no values at abscissae -1
    // and 1.

    return false ;
}


void ParentElement :: Print ()
{
    // Prints the receiver on standard output.

    string s ;

    printf("Parent element %s\n",(GetName(s)).c_str()) ;
}


//______________________________ ParentElementGL ______________________________


Vector<ParentElementGL*>* ParentElementGL :: instances = NULL ;
// Initializes the list of instances (a class variable).


ParentElementGL :: ParentElementGL (int n)
    : ParentElement (n)
{
    // Constructor. Initializes the receiver to a new parent element GL of
    // degree 'n'.
    // Includes the computation of the collocation points and quadrature weights.

    collocationPoints = new RealVector(n+1) ;
    quadratureWeights = new RealVector(n+1) ;
    ZWGJD(collocationPoints,quadratureWeights,ZERO,ZERO) ;
}


ParentElementGL :: ~ParentElementGL ()
{
    // Destructor.
}


void ParentElementGL :: DeleteAllInstances ()
{
    // Deletes the list of the instances of ParentElementGL. This operation can
    // be safely activated only if no object points to any of these parents, 
    // typically if all fields have already been deleted. 
  
    int i ;

    if (instances != NULL) {
	for (i=1 ; i<=instances->GetSize() ; i++)
	    delete instances->At(i) ;
	delete instances ;
	instances = NULL ;
    }
}


FullMatrix* ParentElementGL :: GetDerivMatrix ()
{
    // Returns the derivative matrix of the receiver. Computes this matrix if it
    // does not exist yet.

    int n ;

    if (derivMatrix == NULL) {
	n = polynomialDegree + 1;
	derivMatrix = new FullMatrix(n,n) ;
	DGJD(derivMatrix,this->GetCollocationPoints(),ZERO,ZERO) ;
    }

    return derivMatrix ;
}


Matrix* ParentElementGL :: GetInterpMatrixTo (ParentElement *p)
{
    // Returns the interpolation matrix from the receiver to 'p'.
    // It can be used to calculate u2 = A u1, where u1 contains the discrete
    // values on the receiver of a function u(x), and u2 will contain the 
    // discrete values on 'p' of u(x).

    Matrix     *a ;
    FullMatrix *aa ;  // aa=a, aa defined for avoiding silly compilation troubles

    a = listInterpMatrices->Search(p) ;
    if (a == NULL) {
	if (p == this)
	    a = new IdentityMatrix(GetNbCollocationPoints()) ;
	else {
	    aa = new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints());
	    IGJMD(aa,this->GetCollocationPoints(),p->GetCollocationPoints(),
		  ZERO,ZERO) ;
	    a = aa ;
	}
	listInterpMatrices->Store(p,a) ;
    }

    return a ;
}


string ParentElementGL :: GetName (string s)
{
    // Prints into 's' the name and the polynomial degree of the receiver.
    // Returns 's'.
 
  stringstream buffer;
  buffer << "GL " << polynomialDegree;
  buffer >> s;
  return s ;
}


ParentElement* ParentElementGL :: GetParentElementOfDegree (int n)
{
    // Looks into the list of already created ParentElementGL for the parent of
    // degree 'n'. If such parent does not exist, creates it. 
    // The client should not delete the returned parent.

# ifdef REQUIRE
    Require("acceptable polynomial order 'n'", n >= 0 && n <= MAX_DEGREE) ;
    InMethod("ParentElementGL::GetParentElementOfDegree(n)") ;
# endif

    ParentElementGL *parent ;

    if (instances == NULL) {
	instances = new Vector<ParentElementGL*>(MAX_DEGREE) ;
	instances->SetValues(NULL) ;
    }

    // Note: the parent of degree n is stored at index n+1, in order to allow the
    //       parent of degree 0 to be defined.

    parent = instances->At(n+1) ;
    if (parent == NULL) {
	parent = new ParentElementGL(n) ;
	instances->At(n+1) = parent ;
    }

    return parent ;
}


//____________________________ ParentElementGLL _______________________________


Vector<ParentElementGLL*>* ParentElementGLL :: instances = NULL ;
// Initializes the list of instances (a class variable).


ParentElementGLL :: ParentElementGLL (int n)
    : ParentElement (n)
{
    // Constructor. Initializes the receiver to a GLL parent element of degree
    // 'n'.
    // Includes the computation of the collocation points and quadrature weights.

    collocationPoints = new RealVector(n+1) ;
    quadratureWeights = new RealVector(n+1) ;
    ZWGLJD(collocationPoints,quadratureWeights,ZERO,ZERO) ;
}


ParentElementGLL :: ~ParentElementGLL ()
{
    // Destructor.
}


void ParentElementGLL :: DeleteAllInstances ()
{
    // Deletes the list of the instances of ParentElementGLL. This operation can
    // be safely activated only if no object points to any of these parents, 
    // typically if all fields have already been deleted. 
  
    int i ;

    if (instances != NULL) {
	for (i=1 ; i<=instances->GetSize() ; i++)
	    delete instances->At(i) ;
	delete instances ;
	instances = NULL ;
    }
}


FullMatrix* ParentElementGLL :: GetDerivMatrix ()
{
    // Returns the derivative matrix of the receiver. Computes this matrix if it
    // does not exist yet.

    int n ;

    if (derivMatrix == NULL) {
	n = polynomialDegree + 1;
	derivMatrix = new FullMatrix(n,n) ;
	DGLJD(derivMatrix,this->GetCollocationPoints(),ZERO,ZERO) ;
    }

    return derivMatrix ;
}


Matrix* ParentElementGLL :: GetInterpMatrixTo (ParentElement *p)
{
    // Returns the interpolation matrix from the receiver to 'p'.
    // It can be used to calculate u2 = A u1, where u1 contains the discrete
    // values on the receiver of a function u(x), and u2 will contain the discrete
    // values on 'p' of u(x).

    Matrix     *a ;
    FullMatrix *aa ;  // aa=a, aa defined for avoiding silly compilation troubles

    a = listInterpMatrices->Search(p) ;
    if (a == NULL) {
	if (p == this)
	    a = new IdentityMatrix(GetNbCollocationPoints()) ;
	else {
	    aa =new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints());
	    IGLJMD(aa,this->GetCollocationPoints(),p->GetCollocationPoints(),
		   ZERO, ZERO) ;
	    a = aa ;
	}
	listInterpMatrices->Store(p,a) ;
    }

    return a ;
}


string ParentElementGLL :: GetName (string s)
{
    // Prints into 's' the name and the polynomial degree of the receiver.
    // Returns 's'.

  stringstream buffer;
  buffer << "GLL " << polynomialDegree;
  buffer >> s;
  return s ;
}


ParentElement* ParentElementGLL :: GetParentElementOfDegree (int n)
{
    // Looks into the list of already created ParentElementGLL for the parent of
    // degree 'n'. If such parent does not exist, creates it. 
    // The client should not delete the returned parent.

# ifdef REQUIRE
    Require("acceptable polynomial order 'n'", n >= 1 && n <= MAX_DEGREE) ;
    InMethod("ParentElementGLL::GetParentElementOfDegree(n)") ;
# endif

    ParentElementGLL *parent ;

    if (instances == NULL) {
	instances = new Vector<ParentElementGLL*>(MAX_DEGREE+1) ;
	instances->SetValues(NULL) ;
    }

    parent = instances->At(n) ;
    if (parent == NULL) {
	parent = new ParentElementGLL(n) ;
	instances->At(n) = parent ;
    }

    return parent ;
}


boolean ParentElementGLL :: IsDefinedOnBorders ()
{
    // Returns True, because the receiver has values at abscissae -1 and 1.

    return true ;
}


//__________________________ ParentElementGLLI ________________________________


Vector<ParentElementGLLI*>* ParentElementGLLI :: instances = NULL ;
// Initializes the list of instances (a class variable).


ParentElementGLLI :: ParentElementGLLI (int n)
    : ParentElement (n)
{
    // Constructor. Initializes the receiver to a new parent element GLLI of
    // degree 'n'.
    // Includes the computations of the collocation points.
    // Remember: there is no integration rule (and thus no quadrature weights)
    //           associated with this parent type.

    collocationPoints = new RealVector(n+1) ;
    collocationPointsGLL = new RealVector(n+3) ;  
    quadratureWeights = new RealVector(n+3) ;
    ZWIGLJD(collocationPoints,collocationPointsGLL,quadratureWeights,ZERO,ZERO) ;  
/*
  printf("\n collocationPoints \n");
  collocationPoints->Print();
  collocationPointsGLL->Print();
*/  
}


ParentElementGLLI :: ~ParentElementGLLI ()
{
    // Destructor.
}


void ParentElementGLLI :: DeleteAllInstances ()
{
    // Deletes the list of the instances of ParentElementGLLI. This operation can
    // be safely activated only if no object points to any of these parents, 
    // typically if all fields have already been deleted. 
  
    int i ;

    if (instances != NULL) {
	for (i=1 ; i<=instances->GetSize() ; i++)
	    delete instances->At(i) ;
	delete instances ;
	instances = NULL ;
    }
}


FullMatrix* ParentElementGLLI :: GetDerivMatrix ()
{
/*
// Issues an error message.

NotImplemented("ParentElementGLLI::GetDerivMatrix") ;
return NULL ;
*/

    // Returns the derivative matrix of the receiver. Computes this matrix if it
    // does not exist yet.

    int n ;

    if (derivMatrix == NULL) {
	n = polynomialDegree + 3;
	derivMatrix = new FullMatrix(n,n) ;
	DGLJD(derivMatrix,collocationPointsGLL,ZERO,ZERO) ;
    }

    return derivMatrix ;
  
}


Matrix* ParentElementGLLI :: GetInterpMatrixTo (ParentElement *p)
{
    // Returns the interpolation matrix from the receiver to 'p'.
    // It can be used to calculate u2 = A u1, where u1 contains the discrete
    // values on the receiver of a function u(x), and u2 will contain the discrete
    // values on 'p' of u(x).

    Matrix     *a ;
    FullMatrix *aa ;  // aa=a, aa defined for avoiding silly compilation troubles

//  printf("\n NbCollocationPoints %d \n",this->GetNbCollocationPoints());
    a = listInterpMatrices->Search(p) ;
    if (a == NULL) {
	if (p == this)
	    a = new IdentityMatrix(GetNbCollocationPoints()) ;
	else {
	    aa =new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints());
	    IGLJMD(aa,this->GetCollocationPoints(),p->GetCollocationPoints(),
		   ZERO,ZERO,1) ;
	    a = aa ;
	}
	listInterpMatrices->Store(p,a) ;
    }

    return a ;
}


string ParentElementGLLI :: GetName (string s)
{
    // Prints into 's' the name and the polynomial degree of the receiver.
    // Returns 's'.

  stringstream buffer;
  buffer << "GLLI " << polynomialDegree;
  buffer >> s;
  return s ;
}


ParentElement* ParentElementGLLI :: GetParentElementOfDegree (int n)
{
    // Looks into the list of already created ParentElementGLLI for the parent of
    // degree 'n'. If such parent does not exist, creates it. 
    // The client should not delete the returned parent.

# ifdef REQUIRE
    Require("acceptable polynomial order 'n'", n >= 0 && n <= MAX_DEGREE) ;
    InMethod("ParentElementGLLI::GetParentElementOfDegree(n)") ;
# endif

    ParentElementGLLI *parent ;

    if (instances == NULL) {
	instances = new Vector<ParentElementGLLI*>(MAX_DEGREE) ;
	instances->SetValues(NULL) ;
    }

    // Note: the parent of degree n is stored at index n+1, in order to allow the
    //       parent of degree 0 to be defined.

    parent = instances->At(n+1) ;
    if (parent == NULL) {
	parent = new ParentElementGLLI(n) ;
	instances->At(n+1) = parent ;
    }

    return parent ;
}


//__________________________ ParentElementIGLL ________________________________


Vector<ParentElementIGLL*>* ParentElementIGLL :: instances = NULL ;
// Initializes the list of instances (a class variable).


ParentElementIGLL :: ParentElementIGLL (int n)
    : ParentElement (n)
{
    // Constructor. Initializes the receiver to a new parent element IGLL of
    // degree 'n'.
    // Includes the computations of the collocation points.
    // Remember: there is no integration rule (and thus no quadrature weights)
    //           associated with this parent type.

    collocationPoints = new RealVector(n+1) ;
    ZIGLJD(collocationPoints,ZERO,ZERO) ;
}


ParentElementIGLL :: ~ParentElementIGLL ()
{
    // Destructor.
}


void ParentElementIGLL :: DeleteAllInstances ()
{
    // Deletes the list of the instances of ParentElementIGLL. This operation can
    // be safely activated only if no object points to any of these parents, 
    // typically if all fields have already been deleted. 
  
    int i ;

    if (instances != NULL) {
	for (i=1 ; i<=instances->GetSize() ; i++)
	    delete instances->At(i) ;
	delete instances ;
	instances = NULL ;
    }
}


FullMatrix* ParentElementIGLL :: GetDerivMatrix ()
{
    // Issues an error message.

    NotImplemented("ParentElementIGLL::GetDerivMatrix") ;
    return NULL ;
}


Matrix* ParentElementIGLL :: GetInterpMatrixTo (ParentElement *p)
{
    // Returns the interpolation matrix from the receiver to 'p'.
    // It can be used to calculate u2 = A u1, where u1 contains the discrete
    // values on the receiver of a function u(x), and u2 will contain the discrete
    // values on 'p' of u(x).

    Matrix     *a ;
    FullMatrix *aa ;  // aa=a, aa defined for avoiding silly compilation troubles

    a = listInterpMatrices->Search(p) ;
    if (a == NULL) {
	if (p == this)
	    a = new IdentityMatrix(GetNbCollocationPoints()) ;
	else {
	    aa =new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints());
	    IGLJMD(aa,this->GetCollocationPoints(),p->GetCollocationPoints(),
		   ZERO,ZERO,1) ;
	    a = aa ;
	}
	listInterpMatrices->Store(p,a) ;
    }

    return a ;
}


string ParentElementIGLL :: GetName (string s)
{
    // Prints into 's' the name and the polynomial degree of the receiver.
    // Returns 's'.

  stringstream buffer;
  buffer << "IGLL " << polynomialDegree;
  buffer >> s;
  return s ;
}


ParentElement* ParentElementIGLL :: GetParentElementOfDegree (int n)
{
    // Looks into the list of already created ParentElementIGLL for the parent of
    // degree 'n'. If such parent does not exist, creates it. 
    // The client should not delete the returned parent.

# ifdef REQUIRE
    Require("acceptable polynomial order 'n'", n >= 0 && n <= MAX_DEGREE) ;
    InMethod("ParentElementIGLL::GetParentElementOfDegree(n)") ;
# endif

    ParentElementIGLL *parent ;

    if (instances == NULL) {
	instances = new Vector<ParentElementIGLL*>(MAX_DEGREE) ;
	instances->SetValues(NULL) ;
    }

    // Note: the parent of degree n is stored at index n+1, in order to allow the
    //       parent of degree 0 to be defined.

    parent = instances->At(n+1) ;
    if (parent == NULL) {
	parent = new ParentElementIGLL(n) ;
	instances->At(n+1) = parent ;
    }

    return parent ;
}


RealVector* ParentElementIGLL :: GetWeights ()
{
    // Issues an error message, since this rule has no integration weights
    // (it can be used only as collocation scheme).

    Error("ParentElementIGLL::GetWeights()","no weights (only colloc. points)") ;
    return NULL ;
}


//_____________________________ ParentElementHP _______________________________


Vector<ParentElementHP*>* ParentElementHP :: instances = NULL ;
// Initializes the list of instances (a class variable).


ParentElementHP :: ParentElementHP (int p)
    : ParentElement(p)
{
    // Constructor. Creates a h-p parent element of degree 'p'.
}


ParentElementHP :: ~ParentElementHP ()
{
    // Destructor.
}


void ParentElementHP :: DeleteAllInstances ()
{
    // Deletes the list of the instances of ParentElementHP. This operation can
    // be safely activated only if no object points to any of these parents, 
    // typically if all fields have already been deleted. 
  
    int i ;

    if (instances != NULL) {
	for (i=1 ; i<=instances->GetSize() ; i++)
	    delete instances->At(i) ;
	delete instances ;
	instances = NULL ;
    }
}


real ParentElementHP :: EvalDerivHP (int j, real r)
{
    // Returns the derivative of the 'j'-th shape function of the receiver,
    // evaluated at 'r'. Normally, 'r' is in [-1.,1].

# ifdef REQUIRE
    Require("valid argument 'j'", j >= 0 && j <= polynomialDegree) ;
    InMethod("ParentElementHP::EvalDerivHP(j,r)") ;
# endif

    if (j == 0)
	return -HALF ;

    else if (j == polynomialDegree)
	return HALF ;

    else 
	return sqrt(HALF*(j+j+ONE)) * EvalLegendre(j,r) ;
}


real ParentElementHP :: EvalHP (int j, real r)
{
    // Returns the 'j'-th shape function of the receiver evaluated at 'r'.
    // Normally, 'r' is in [-1.,1].

# ifdef REQUIRE
    Require("valid argument 'j'", j >= 0 && j <= polynomialDegree) ;
    InMethod("ParentElementHP::EvalHP(j,r)") ;
# endif

    real factor, pPlus, pMinus ;

    if (j == 0)
	return HALF*(ONE-r) ;

    else if (j == polynomialDegree)
	return HALF*(ONE+r) ;

    else {
	factor = ONE / sqrt(FOUR*j+TWO) ;
	pPlus  = EvalLegendre(j+1,r) ;                   // P_j+1(r)
	pMinus = EvalLegendre(j-1,r) ;                   // P_j-1(r)
	return factor*(pPlus-pMinus) ;
    }
}


FullMatrix* ParentElementHP :: GetDerivMatrix ()
{
    // Issues an error message, because h-p formulations have no collocation
    // points.
    // You must use instead method GetInterpDerivMatrixTo(p).

    NotImplemented("ParentElementHP::GetDerivMatrix()") ;
    return NULL ;
}


FullMatrix* ParentElementHP :: GetInterpDerivMatrixTo (ParentElement* p)
{
    // Returns the interpolated derivative matrix from the receiver to 'p'.
    // Every coefficient D(i,j) contains the derivatives of the (j+1)-th shape
    // function of the receiver at the (i+1)-th collocation point of 'p'.

    FullMatrix *a ;
    RealVector *z ;
    int        i, j ;

    a = (FullMatrix*) listInterpDerivMatrices->Search(p) ;
    if (a == NULL) {
	a = new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints()) ;
	z = p->GetCollocationPoints() ;
	for (j=0 ; j<=polynomialDegree ; j++)
	    for (i=0 ; i<=p->GetPolynomialDegree() ; i++)
		a->At(i+1,j+1) = EvalDerivHP(j,z->At(i+1)) ;
	listInterpDerivMatrices->Store(p,a) ;
    }

    return a ;
}


Matrix* ParentElementHP :: GetInterpMatrixTo (ParentElement* p)
{
    // Returns the interpolation matrix from the receiver to 'p'.
    // It can be used to calculate u2 = A u1, where u1 contains the discrete
    // values on the receiver of a function u(x), and u2 will contain the discrete
    // values on 'p' of u(x).

# ifdef REQUIRE
    Require("valid argument 'p'", p->GetType() != HP_PT) ;
    InMethod("ParentElementHP::GetInterpMatrixTo(p)") ;
# endif

    Matrix     *a ;
    FullMatrix *aa ;
    RealVector *z ;
    int        i, j ;

    a = listInterpMatrices->Search(p) ;
    if (a == NULL) {
	if (p == this)
	    a = new IdentityMatrix(GetNbCollocationPoints()) ;
	else {
	    aa=new FullMatrix(p->GetNbCollocationPoints(),GetNbCollocationPoints());
	    z = p->GetCollocationPoints() ;
	    for (j=0 ; j<=polynomialDegree ; j++)
		for (i=0 ; i<=p->GetPolynomialDegree() ; i++)
		    aa->At(i+1,j+1) = EvalHP(j,z->At(i+1)) ;
	    a = aa ;
	}
	listInterpMatrices->Store(p,a) ;
    }

    return a ;
}


string ParentElementHP :: GetName (string s)
{
    // Prints into 's' the name and the polynomial degree of the receiver.
    // Returns 's'.

  stringstream buffer;
  buffer << "HP " << polynomialDegree;
  buffer >> s;
  return s ;
}


ParentElement* ParentElementHP :: GetParentElementOfDegree (int n)
{
    // Looks into the list of already created ParentElementHP for the parent of
    // degree 'n'. If such parent does not exist, creates it. 
    // The client should not delete the returned parent.

# ifdef REQUIRE
    Require("acceptable polynomial order 'n'", n >= 1 && n <= MAX_DEGREE) ;
    InMethod("ParentElementHP::GetParentElementOfDegree(n)") ;
# endif

    ParentElementHP *parent ;

    if (instances == NULL) {
	instances = new Vector<ParentElementHP*>(MAX_DEGREE+1) ;
	instances->SetValues(NULL) ;
    }

    parent = instances->At(n) ;
    if (parent == NULL) {
	parent = new ParentElementHP(n) ;
	instances->At(n) = parent ;
    }

    return parent ;
}


boolean ParentElementHP :: IsDefinedOnBorders ()
{
    // Returns True, because the receiver has values at abscissae -1 and 1.

    return true ;
}


