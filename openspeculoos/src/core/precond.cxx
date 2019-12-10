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

// precond.cxx

#include "core/precond.hxx"
#include "core/problem.hxx"
#include "core/mesh.hxx"
#include "core/field.hxx"
#include "core/element.hxx"
#include "core/parent.hxx"
#include "core/matrix.hxx"
#include "core/timer.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <math.h>


//________________________________ Preconditioner _____________________________


Preconditioner :: Preconditioner () 
    : GarbagedObject ()
{
    // Constructor.
  
    problem = NULL ;
}


Preconditioner :: ~Preconditioner ()
{
    // Destructor.
}


void Preconditioner :: SetProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the problem to which the receiver is attached. 'prob' is
    // used to give the operator A that the receiver will precondition.
    // 'prob' can be NULL.

    problem = prob ;
}


//____________________________ FieldPreconditioner ____________________________


FieldPreconditioner :: FieldPreconditioner ()
    : Preconditioner ()
{
    // Constructor.

    inverseP = NULL ;
}


FieldPreconditioner :: ~FieldPreconditioner ()
{
    // Destructor.

    unref(inverseP) ;
}


void FieldPreconditioner :: Precondition (Field* r, Field* z)
{
    // Solves for z the system P z = r. 

# ifdef REQUIRE
    Require("z has same interpolation as r", z->HasSameInterpolationAs(r)) ;
    Require("z has same interpolation as u", problem->GetUnknown()
	    ->HasSameInterpolationAs(r));
    Require("z has same nb components as r", z->GetNbComponents() ==
	    r->GetNbComponents()) ;
    InMethod("FieldPreconditioner::Precondition(r,z)") ;
# endif

    z->CopyFrom(r) ;
    z->Multiply(GetInverseP()) ;
}


void FieldPreconditioner :: SetProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the problem to which the receiver is attached.
    // 'prob' may be NULL

    problem = prob ;
    unref(inverseP) ;
}


//___________________________ DiagonalPreconditioner __________________________


Field* DiagonalPreconditioner :: GetInverseP ()
{
    // Returns inv(P), where P is the diagonal of A.
    // For every i such that P(i,i)=0 (typically, if the i-th dof is subject to a
    // Dirichlet condition), sets inv(P)(i,i) to 1.
    // Programming note. As opposed to classes EuclidianNormPreconditioner and
    //    InfiniteNormPreconditioner, and in order to speed up the calculations,
    //    here the loop cycles on the elements, not on every dof.

    Field     *u ;
    FlatField *temp1, *temp2 ;
    int       i, maxDof ;

    if (! inverseP) {

	Printf("\nStart computing inverse diagonal operator ...\n") ;

	u        = problem->GetUnknown() ;  
	inverseP = u->FieldDuplicateEmpty() ;
	temp1    = u->GetMain()->DuplicateEmpty() ;
	temp2    = u->GetMain()->DuplicateEmpty() ;
	maxDof   = u->GetMain()->GetMaxNbLocalDofs() ;

	for (i=1 ; i<=maxDof ; i++) {
	    if (i>1)
		temp1->SetLocalDofs(i-1,ZERO) ;
	    temp1->SetLocalDofs(i,ONE) ;

	    problem->GetFlatAx(temp1,temp2) ;

	    inverseP->GetMain()->SetLocalDofs(i,temp2) ;
	}

	delete temp1 ;
	delete temp2 ;

	inverseP->Assemble() ;
	inverseP->Inverse() ;

	Printf("OK inverse diagonal operator\n\n") ;
    }

    return inverseP ;
}


//_____________________________ MassPreconditioner ____________________________


Field* MassPreconditioner :: GetInverseP ()
{
    // Returns inv(P), where P is the assembled mass matrix of the problem.
    // inv(P) has one component.

    if (! inverseP) {
	inverseP = problem->GetInverseB() ;
	ref(inverseP) ;
    }

    return inverseP ;
}


void MassPreconditioner :: Precondition (Field* r, Field* z)
{
    // Solves for z the diagonal system P z = r. 

# ifdef REQUIRE
    Require("z has same interpolation as r", z->HasSameInterpolationAs(r)) ;
    Require("z has same interpolation as u", problem->GetUnknown()
	    ->HasSameInterpolationAs(r));
    Require("z has same nb components as r", z->GetNbComponents() ==
	    r->GetNbComponents()) ;
    InMethod("MassPreconditioner::Precondition(r,z)") ;
# endif

    int i ;

    z->CopyFrom(r) ;

    GetInverseP () ;
    for (i=1 ; i <= z->GetNbComponents() ; i++)
	z->Multiply(inverseP,1,i) ;
}


//_________________________ EuclidianNormPreconditioner _______________________


Field* EuclidianNormPreconditioner :: GetInverseP ()
{
    // Returns inv(P), where P is the diagonal matrix made of the euclidian norm
    // of the i-th column of A.
    // For every i such that P(i,i)=0 (typically, if the i-th dof is subject to a
    // Dirichlet condition), sets inv(P)(i,i) to 1.
    // Since inv(P) is diagonal, it can be stored as a field.

    Field *u, *temp1, *temp2 ;
    real  norm ;
    int   i, nbdof ;

    if (! inverseP) {

	printf("\nStart computing Euclidian-norm preconditioner ...\n") ;

	u        = problem->GetUnknown() ;
	inverseP = u->FieldDuplicateEmpty() ;
	temp1    = u->FieldDuplicateEmpty() ;
	temp2    = u->FieldDuplicateEmpty() ;

	nbdof    = u->GetNbDofs() ;

	for (i=0 ; i<nbdof ; i++) {

	    if (i>0)
		temp1->SetDof(i-1,ZERO) ;             // set all dofs to 0,
	    temp1->SetDof(i,ONE) ;                  // except the i-th dof (set to 1)

	    problem->GetAx(temp1,temp2) ;           // temp2 = i-th column of A

	    norm = temp2->GetEuclidianNorm() ;
	    if (fabs(norm) < 1.e-10)
		norm = ONE ;
	    else  
		norm = ONE / norm ;
	    inverseP->SetDof(i,norm) ;
	}

	printf("done\n") ;
	delete temp1 ;
	delete temp2 ;
    }

    return inverseP ;
}


//_________________________ InfiniteNormPreconditioner ________________________


Field* InfiniteNormPreconditioner :: GetInverseP ()
{
    // Returns inv(P), where P is the diagonal matrix made of the infinite norm
    // of the i-th column of A.
    // For every i such that P(i,i)=0 (typically, if the i-th dof is subject to a
    // Dirichlet condition), sets inv(P)(i,i) to 1.
    // Since inv(P) is diagonal, it can be stored as a field.

    Field *u, *temp1, *temp2 ;
    real  norm ;
    int   i, nbdof ;

    if (! inverseP) {

	printf("\nStart computing infinite-norm preconditioner ...\n") ;

	u        = problem->GetUnknown() ;
	inverseP = u->FieldDuplicateEmpty() ;
	temp1    = u->FieldDuplicateEmpty() ;
	temp2    = u->FieldDuplicateEmpty() ;

	nbdof    = u->GetNbDofs() ;

	for (i=0 ; i<nbdof ; i++) {

	    if (i>0)
		temp1->SetDof(i-1,ZERO) ;             // set all dofs to 0,
	    temp1->SetDof(i,ONE) ;                  // except the i-th dof (set to 1)

	    problem->GetAx(temp1,temp2) ;           // temp2 = i-th column of A

	    norm = temp2->GetInfiniteNorm() ;
	    if (fabs(norm) < 1.e-10)
		norm = ONE ;
	    else  
		norm = ONE / norm ;
	    inverseP->SetDof(i,norm) ;
	}

	printf("done\n") ;
	delete temp1 ;
	delete temp2 ;
    }

    return inverseP ;
}


//___________________________ PseudoLaplacePreconditioner _____________________


PseudoLaplacePreconditioner :: PseudoLaplacePreconditioner (CouzyType typ)
    : Preconditioner()
{
    // Constructor.

    type               = typ ;
    elementaryMatrices = NULL ;
}


PseudoLaplacePreconditioner :: ~PseudoLaplacePreconditioner ()
{
    // Destructor.

    int i ;

    if (elementaryMatrices)
	for (i=1 ; i<=elementaryMatrices->GetSize() ; i++)
	    delete elementaryMatrices->At(i) ;

    delete elementaryMatrices ;
}


void PseudoLaplacePreconditioner :: Precondition (Field* r, Field* z)
{
    // Solves for z the P z = r, i.e., initializes z to
    //   z = inv (D inv(B) D(T)) r,
    // computed element by element, i.e., approximately.
    //
    // On every element, the boundary conditions on B (the mass matrix on the
    // velocity grid) are:
    // - if type = DIRICHLET, homogeneous Dirichlet conditions,
    // - if type = DIRICHLET_NEUMANN, homogeneous Dirichlet conditions on the 
    //   part of the element that is subject to Dirichlet conditions, homogeneous
    //   Neumann conditions on the rest (i.e., the interface with other 
    //   elements).
    //
    // r is not modified.

# ifdef REQUIRE
    Require("same type (r,unknown)", r->HasSameTypeAs(problem->GetUnknown())) ;
    Require("same type (z,unknown)", z->HasSameTypeAs(problem->GetUnknown())) ;
    Require("r has same interpolation as p", r->HasSameInterpolationAs(z)) ;
    Require("z has same interpolation as p", z->HasSameInterpolationAs(
		problem->GetUnknown())) ;
    Require("z has same nb components as r", z->GetNbComponents() ==
	    r->GetNbComponents()) ;
    InMethod("PseudoLaplacePreconditioner::Precondition(r,z)") ;
# endif

    Mesh         *mesh ;
    Element      *elem ;
    FullMatrixLU *E ;
    RealVector   *rr, *zz ;
    int          k;

    // 1. Get the ingredients of inv(E)

    if (! elementaryMatrices)
	SetElementaryMatrices() ;

    // 2. Initialize z to inv(E).r

    mesh = problem->GetMesh() ;

    std::vector<int> const& localElements = mesh->GetLocalElements();
    for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
        k = localElements[iLocal];
	elem = mesh->At(k) ;
        E  = elementaryMatrices->At(k) ;
        rr = elem->GetElementaryField(r->GetMain()->GetIndex())->GetComponent(1);
        zz = elem->GetElementaryField(z->GetMain()->GetIndex())->GetComponent(1);
        zz->CopyFrom(rr) ;
        E->Solve(zz) ;
    }
}


void PseudoLaplacePreconditioner :: SetElementaryMatrices ()
{
    // Computes the elementary matrices E = D inv(B) D(T).

# ifdef REQUIRE
    Require("problem is steady Stokes", problem->IsSteadyStokes()) ;
    InMethod("PseudoLaplacePreconditioner::SetElementaryMatrices()") ;
# endif

    SteadyStokes  *prob ;
    Mesh          *mesh ;
    MortaredField *temp ;
    FlatField     *inverseB, *p, *v, *rule, *homog ;
    Element       *elem ;
    FullMatrixLU  *E ;
    RealVector    *valP ;
    int           i, j, k, nbDofs, maxNbDofs, nbComp ;
    real          norm ;

    Printf("\nStart computing pseudo-Laplace preconditioner ...\n") ;

    // 1. Get inv(B) with adequate boundary conditions. B is not assembled.

    prob     = (SteadyStokes*)problem ;
    inverseB = prob->GetVelocity()->GetMain()->DuplicateEmpty() ;
    nbComp   = inverseB->GetNbComponents() ;
    for (i=1 ; i<=nbComp ; i++)
	inverseB->CopyInterpolateFrom(prob->GetMesh()->GetJacobian(),1,i) ;
    inverseB->MultiplyByWeights() ;
    inverseB->Inverse() ;
    norm = inverseB->GetEuclidianNorm() ;
    if (Mpi_rank==0) printf(" step1 : norm inverseB = %lf \n",norm) ;
    if (type == DIRICHLET) {   
	// homogenous Dirichlet on the whole contour of every element

	temp  = prob->GetVelocity()->DuplicateEmpty() ;
	homog = temp->GetMortar()->GetMain() ;     // defined on the whole skeleton
	inverseB->CopyValuesFrom(homog,POINTWISE_NC) ; 
	delete temp ;
    }
    else {  // type == DIRICHLET_NEUMANN
	// homogeneous Dirichlet on the external contour of every element,
	// homogeneous Neumann (i.e., nothing) on the interface of every element
	prob->ApplyHomogeneousDirichletConditionsV(inverseB) ;
    }
    norm = inverseB->GetEuclidianNorm() ;
    if (Mpi_rank==0) printf(" Precond SetEltMatrice b4 phase2 : norm inverseB = %lf \n", norm) ;
    // 2. Compute E = D inv(B) D(T) column by column, by probing.
    //    This is done simultaneously for all elements.

    p    = prob->GetPressure()->DuplicateEmpty() ;
    v    = prob->GetVelocity()->GetMain()->DuplicateEmpty() ;
    rule = prob->GetIntegrationRuleMass() ;

    mesh               = prob->GetMesh() ;
    elementaryMatrices = new Vector<FullMatrixLU*>(mesh->GetSize()) ;
    maxNbDofs          = p->GetMaxNbLocalDofs() ;

    for (j=1 ; j<=maxNbDofs ; j++) {         // j-th column of E of every element

	// set all dofs to 0, except the j-th dof of every element

	p->SetValues(ZERO) ;
	p->SetLocalDofs(j,ONE) ;

	// compute E p

	v->SetToWeakGradient(p,rule) ;               // - D(T) p
	v->Multiply(inverseB) ;                      // - inv(B) D(T) p 
	p->SetToWeakDivergence(v,rule) ;             // - D inv(B) D(T) p = - E p

        std::vector<int> const& localElements = mesh->GetLocalElements();
        for (int iLocal=0; iLocal<localElements.size(); ++iLocal) {
            k = localElements[iLocal];
            elem = mesh->At(k) ;
            valP   = elem->GetElementaryField(p->GetIndex())->GetComponent(1) ;
            nbDofs = valP->GetSize() ;

            if (j <= nbDofs) {       // not all elements have maxNbDofs dofs

                if (j == 1) {                              // the first column
                    E = new FullMatrixLU(nbDofs,nbDofs) ;
                    elementaryMatrices->At(k) = E ;
                }
                else                                       // a subsequent column
                    E = elementaryMatrices->At(k) ;

                for (i=1 ; i<=nbDofs ; i++)
                    E->At(i,j) = - valP->At(i) ;
            }
	}
    }

    delete inverseB ;
    delete p ;
    delete v ;

    Printf("OK pseudo-Laplace preconditioner\n\n") ;
}


//____________________________ CouzyPreconditioner ____________________________


CouzyPreconditioner :: CouzyPreconditioner (CouzyType typ)
    : Preconditioner()
{
    // Constructor. Initializes the receiver to a Couzy preconditioner with
    // homogeneous boundary conditions on the velocity mass operator of type
    // 'typ' (Neumann or Dirichlet).

    pseudoLaplace = new PseudoLaplacePreconditioner(typ) ;
}


CouzyPreconditioner :: ~CouzyPreconditioner ()
{
    // Destructor.

    delete pseudoLaplace ;
}


void CouzyPreconditioner :: Precondition (Field* r, Field* z)
{
    // Solves for z the diagonal system P z = r, i.e., initializes z to
    //   z = {nu inv(B~) + lambda inv(E)} r. 
    //
    // r is not modified.

# ifdef REQUIRE
    Require("same type (r,unknown)", r->HasSameTypeAs(problem->GetUnknown())) ;
    Require("same type (z,unknown)", z->HasSameTypeAs(problem->GetUnknown())) ;
    Require("r has same interpolation as p", r->HasSameInterpolationAs(z)) ;
    Require("z has same interpolation as p", z->HasSameInterpolationAs(
		problem->GetUnknown())) ;
    Require("z has same nb components as r", z->GetNbComponents() ==
	    r->GetNbComponents()) ;
    InMethod("CouzyPreconditioner::Precondition(r,z)") ;
# endif

    Field *z1 ;
    real  c1, c2 ;

    // 1. Get c1 and c2

    c1 = problem->GetC1() ;
    c2 = problem->GetC2() ;

    // 2. Compute z1 = inv(B) r

    if (c1 != ZERO) {
	z1 = z->FieldGetWork() ;
	z1->CopyFrom(r) ;
	problem->ApplyInverseB(z1) ;
    }

    // 3. Compute z2 = inv(E) r

    pseudoLaplace->Precondition(r,z) ;     // usually c2 is not 0

    // 4. Compute z = c1 z1 + c2 z2

    if (c1 != ZERO) {                      // if c1 = 0, z need not be multiplied
	z->MultiplyAndAdd(c2,z1,c1) ;        //                               by c2
	z->Retrieve(z1) ;
    }
}


void CouzyPreconditioner :: SetProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the problem of the receiver.

# ifdef REQUIRE
    Require("'prob' is a SteadyStokes problem", prob->IsSteadyStokes()) ;
    InMethod("CouzyPreconditioner::SetProblem(prob)") ;
# endif

    problem = prob ;

    pseudoLaplace->SetProblem(prob) ;
}


