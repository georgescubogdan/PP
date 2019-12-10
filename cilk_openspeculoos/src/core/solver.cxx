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

// solver.cxx

#include "core/solver.hxx"
#include "core/problem.hxx"
#include "core/precond.hxx"
#include "core/field.hxx"
#include "core/util.hxx"
#include "core/timer.hxx"
#include "core/parallel.hxx"
#include "core/options.hxx"
#include <limits>
#include <math.h>


// Ajout Vincent Keller (needs libmonitor)
//#include "MPIMonitor_v1.1/mpimonitor.h"


enum {CONVERGENCE, DIVERGENCE, NO_CONVERGENCE, STAGNATION} ;


//________________________________ Solver _____________________________________


Solver :: Solver ()
    : GarbagedObject ()
{
    // Constructor.

    name       = "(no name)" ;
    verbosity  = 0 ;
}

Solver :: ~Solver ()
{
    // Destructor.

}

void Solver :: Print ()
{
    // Default implementation: does nothing.
}

void Solver :: PrintTime ()
{
    // Default implementation: does nothing.
}

void Solver :: PrintTime (FILE *file)
{
    // Default implementation: does nothing.
}

void Solver :: SetAdjustTolerance (boolean)
{
    // Default implementation: does nothing.
}

void Solver :: SetName (string s)
{
    // Sets the receiver's name to 's'. Used for printings.

    name = s ;
}

void Solver :: SetVerbosity (int v)
{
    // Sets the verbosity of the receiver to 'v' (0 for mute).

    verbosity = v ;
}


//____________________________ ConjugateGradient ______________________________


ConjugateGradient :: ConjugateGradient (int maxIter, real tol)
    : Solver ()
{
    // Constructor. Initializes the receiver to a conjugate-gradient solver, with
    // at most 'maxIter' iterations and 'tol' as tolerance on the residue.

    maxIterations   = maxIter ; 
    tolerance       = tol ;
    normRHS         = UNINITIALIZED_REAL ;
    adjustTolerance = false ;
    iter            = 0 ;
    iterCumulated   = 0 ;
}

ConjugateGradient :: ~ConjugateGradient ()
{
    // Destructor.
}

void ConjugateGradient :: Print ()
{
    // Prints on standard output the number of performed iterations.

  
    Printf("\nCG '%s': nb iterations: %d  cumulated: %d\n",
	   name.c_str(),iter,iterCumulated) ;
}

void ConjugateGradient :: SetAdjustTolerance (boolean b)
{
    // If 'adjustTolerance' is 'true', then the convergence criterion will be
    // relaxed (see header file). Default value: 'false'.

    adjustTolerance = b ;
    normRHS         = UNINITIALIZED_REAL ;
} 

void ConjugateGradient :: Solve (SemiDiscreteProblem* prob)
{
    // Solves 'prob', i.e., initializes the unknown field "u" of 'prob' to the
    // solution of the system A u = RHS of 'prob'.

    Field          *unknown, *p, *q, *r, *z ;
    Preconditioner *precond ;
    real           normR, ratio, denom, alpha, beta, rho, rho1=ONE ;
    int            result ;



    // 1. Compute the right-hand side 'r' and its norm (used as a reference norm 
    //    for the iterative algorithm)

    unknown = prob->GetUnknown() ;
    r       = unknown->FieldGetWork() ;

    prob->GetRHS(r) ;

    if (adjustTolerance == false || normRHS == UNINITIALIZED_REAL)
	normRHS = r->GetEuclidianNorm() ;

    if (verbosity >= 3) {
	Printf(" iter%4d  ratio:%14.7e  norm(RHS)=%14.7e\n", 0, 1., normRHS) ;
    }

    // 2. Conjugate-gradient iterations

    precond = prob->GetPreconditioner() ;
    p       = unknown->FieldGetWork() ;
    q       = unknown->FieldGetWork() ;
    z       = q ;
    unknown->SetValues(ZERO) ;


    // Code introduced by Jonas Latt: This condition is used to
    //   avoid solving a singular system when normRHS is too close
    //   to zero.
    if (normRHS <= std::numeric_limits<double>::epsilon() ) {
	iter   = 0 ;
	result = CONVERGENCE ;
    }

    else {
	result = NO_CONVERGENCE ;

	for (iter=1; iter<=maxIterations ; iter++) {
	    iterCumulated++ ;

	    if (precond)
		precond->Precondition(r,z) ;
	    else
		z->CopyFrom(r) ;

	    rho = r->Dot(z) ;

	    if (iter == 1)
		p->CopyFrom(z) ;
	    else {
		beta = rho / rho1 ;
		p->MultiplyAndAdd(beta,z) ;
	    }

	    prob->GetAx(p,q) ;

	    denom = p->Dot(q) ;

	    if (fabs(denom) < 1.e-30 * fabs(rho)) {
	      Printf("\nCG '%s': denom = %.4e  rho = %.4e\n",name.c_str(),denom,rho) ;
		Error("ConjugateGradient::Solve","indefinite system") ;
	    }

	    alpha = rho / denom ;
	    unknown->Add(p,alpha) ;
	    r->Add(q,-alpha) ;

	    normR = r->GetEuclidianNorm() ;
	    ratio = normR / normRHS ;
	    if (verbosity >= 3 && iter%100==0) {
		Printf(" iter%4d  ratio:%14.7e\n", iter, ratio) ;
	    }
	    if (verbosity == 3)
		Printf(" iter%4d  ratio:%14.7e     normR =%14.7e\n", iter, ratio, normR) ;
      

	    if (ratio <= tolerance) {
		result = CONVERGENCE ;
		break ;
	    }
	    else if (ratio > VLARGE) {
		result = DIVERGENCE ;
		break ;
	    }
 
	    rho1 = rho ;
	}
    }

    // 3. Print convergence information

    if (result == CONVERGENCE) {
	if (verbosity >= 2)
	    Printf(" iter%4d  ratio:%14.7e   normR =%14.7e\n", iter, ratio, normR) ;
	if (verbosity >= 1) {
	  Printf("Iterative solver '%s' converged in %d iterations\n",name.c_str(),iter) ;
	}
    }
    else if (result == DIVERGENCE) {
      Printf("Iterative solver '%s' diverged at iteration %d\n",name.c_str(),iter) ;
    }
    else {
	Printf("Iterative solver '%s' not converged after %d iterations\n",
	       name.c_str(),maxIterations) ;
    }

    // 4. Release the work fields

    unknown->Retrieve(p) ;
    unknown->Retrieve(q) ;
    unknown->Retrieve(r) ;


}

//____________________________ GMRES ______________________________


GMRES :: GMRES (int maxIter, real tol)
    : Solver ()
{
    // Constructor. Initializes the receiver to a conjugate-gradient solver, with
    // at most 'maxIter' iterations and 'tol' as tolerance on the residue.

    maxIterations   = maxIter ; 
    tolerance       = tol ;
    normRHS         = UNINITIALIZED_REAL ;
    adjustTolerance = false ;
    iter            = 0 ;
    iterCumulated   = 0 ;
}

GMRES :: ~GMRES ()
{
    // Destructor.
}

void GMRES :: Print ()
{
    // Prints on standard output the number of performed iterations.

  
    Printf("\nCG '%s': nb iterations: %d  cumulated: %d\n",
	   name.c_str(),iter,iterCumulated) ;
}

void GMRES :: SetAdjustTolerance (boolean b)
{
    // If 'adjustTolerance' is 'true', then the convergence criterion will be
    // relaxed (see header file). Default value: 'false'.

    adjustTolerance = b ;
    normRHS         = UNINITIALIZED_REAL ;
} 

void GMRES :: Solve (SemiDiscreteProblem* prob)
{
// subroutine gmres(x,b,restart,tol,maxit,n,flag,relres)

//   use Constants

//   implicit none

//   logical restarted
//   integer maxit,n,restart,inner,outer,imin,jmin,i,j,k,l,flag,ipiv(n)
//   double precision b(n),x(n),x0(n),xmin(n),r(n),iter(2),vh(n),vh1(n),u(n),u1(n),u2(n),vrf(n)
//   double precision tol,relres,n2b,tolb,normrmin,rho,phibar,rt,temp,c,s,normr,norm

//   double precision, dimension(:),   allocatable :: resvec
//   double precision, dimension(:,:), allocatable :: V
//   double precision, dimension(:),   allocatable :: h
//   double precision, dimension(:,:), allocatable :: QT
//   double precision, dimension(:,:), allocatable :: RM
//   double precision, dimension(:),   allocatable :: f
//   double precision, dimension(:,:), allocatable :: W
//   double precision, dimension(:),   allocatable :: q
//   double precision, dimension(:),   allocatable :: y
//   double precision, dimension(:),   allocatable :: tempv


//   if (restart .eq. n) then
//      restarted = zero
//   else
//      restarted = one
//   end if

//   if (restarted) then
//      outer = maxit
//      inner = restart
//   else
//      outer = one
//      inner = maxit
//   end if

//   // Check for all zero right hand side vector => all zero solution

//   n2b = dsqrt(dot_product(b,b))            // Norm of rhs vector, b
//   if (n2b .eq. zero) then                  // if    rhs vector is all zeros
//      x = dzero                             // then  solution is all zeros
//      flag = one                            // a valid solution has been obtained
//      relres = dzero                        // the relative residual is actually 0/0
//      iter(1) = zero ; iter(2) = zero       // no iterations need be performed
//      allocate(resvec(1))
//      resvec = dzero                        // resvec(1) = norm(b-A*x) = norm(0)
//      return
//   end if

//   x0 = x


//   // Set up for the method
//   flag = one
//   xmin = x                          // Iterate which has minimal residual so far
//   imin = zero                       // "Outer" iteration at which xmin was computed
//   jmin = zero                       // "Inner" iteration at which xmin was computed
//   tolb = tol * n2b                  // Relative tolerance


//   call GetH(x,r) ;  r = b-r         // r = b - matmul(A,x)

//   normr = dsqrt(dot_product(r,r))   // Norm of residual

//   if (normr .le. tolb)  then        // Initial guess is a good enough solution
//      flag = zero
//      relres = normr / n2b
//      iter(1) = zero ; iter(2) = zero
//      allocate(resvec(1))
//      resvec = normr
//      deallocate(resvec)
//      return
//   end if

//   allocate(resvec(inner*outer+1))

//   resvec    = dzero                 // Preallocate vector for norm of residuals
//   resvec(1) = normr                 // resvec(1) = norm(b-A*x0)
//   normrmin  = normr                 // Norm of residual from xmin
//   rho       = done


//   // loop over "outer" iterations (unless convergence or failure)

//   allocate(V(n,inner+1))        ; V  = dzero    // Arnoldi vectors
//   allocate(h(inner+1))          ; h  = dzero    // upper Hessenberg st A*V = V*H ...
//   allocate(QT(inner+1,inner+1)) ; QT = dzero    // orthogonal factor st QT*H = RM
//   allocate(RM(inner,inner))     ; RM = dzero    // upper triangular factor st H = Q*RM
//   allocate(f(inner))            ; f  = dzero    // y = RM\f => x = x0 + V*y
//   allocate(W(n,inner))          ; W  = dzero    // W = V*inv(RM)
//   allocate(q(inner))            ; q  = zero     // "inner" iteration counter
//   allocate(y(inner))
//   allocate(tempv(inner))


//   out : do i = 1, outer
     
//      vh1 = r ;  // vh1 = M1 \ r
//      vh = vh1 ; // vh = M2 \ vh1
//      h(1) = dsqrt(dot_product(vh,vh))
//      V(:,1) = vh / h(1)
//      QT(1,1) = done
//      phibar = h(1)

//      in : do j = 1, inner
//         call GetH(V(:,j),u2)    //u2 = matmul(H,V(:,j))
//         u1 = u2 ;   // u1 = M1 \ u2
//         u = u1 ;    // u = M2 \ u1

//         do k = 1, j
//            h(k) = dot_product(V(:,k),u)
//            u = u - h(k) * V(:,k)
//         end do
//         h(j+1) = dsqrt(dot_product(u,u))
//         V(:,j+1) = u / h(j+1)
//         RM(1:j,j) = matmul(QT(1:j,1:j),h(1:j))
//         rt = RM(j,j)

//         // find cos(theta) and sin(theta) of Givens rotation
//         if (h(j+1) .eq. dzero) then
//            c = done                      // theta = 0
//            s = dzero
//         else if (abs(h(j+1)) .gt. abs(rt)) then
//            temp = rt / h(j+1)
//            s = 1.d0 / dsqrt(1.d0 + abs(temp)**2.d0) // pi/4 < theta < 3pi/4
//            c = - temp * s
//         else
//            temp = h(j+1) / rt
//            c = 1.d0 / dsqrt(1.d0 + abs(temp)**2.d0) // -pi/4 <= theta < 0 < theta <= pi/4
//            s = - temp * c
//         end if
//         RM(j,j) = c * rt - s * h(j+1)

//         q(1:j) = QT(j,1:j)
//         QT(j,1:j) = c*q(1:j)
//         QT(j+1,1:j) = s*q(1:j)
//         QT(j,j+1) = -s
//         QT(j+1,j+1) = c
//         f(j) = c * phibar
//         phibar = s * phibar

//         if (j .lt. inner) then

//            W(:,j) = (V(:,j)-matmul(W(:,1:j-1),RM(1:j-1,j)))/RM(j,j)

//            x = x + f(j) * W(:,j)                   // form the new inner iterate
//         else // j == inner
//            tempv(i:j) = f(1:j)
//            call dgesv(j,one,RM(1:j,1:j),j,ipiv(1:j),tempv(1:j),j,flag)
//            vrf = matmul(V(:,1:j),tempv(1:j))
//            x = x0 + vrf                            // form the new outer iterate
//         end if

//         call GetH(x,r) ; r = b-r
//         normr = dsqrt(dot_product(r,r))
//         resvec((i-1)*inner+j+1) = normr
//         if (normr .le. tolb) then                  // check for convergence
//            if (j .lt. inner) then
//               y(1:j) = f(1:j)
//               call dgesv(j,one,RM(1:j,1:j),j,ipiv(1:j),y(1:j),j,flag)
//               x = x0 + matmul(V(:,1:j),y(1:j))     // more accurate computation of xj
//               call GetH(x,r) ; r = b-r             // r = b - matmul(H,x)
//               normr = dsqrt(dot_product(r,r))
//               resvec((i-1)*inner+j+1) = normr
//            end if
//            if (normr .le. tolb) then               // check using more accurate xj
//               flag = zero
//               iter(1) = i ; iter(2) = j
//               exit in
//            end if
//         end if

//         if (normr .lt. normrmin) then              // update minimal norm quantities
//            normrmin = normr
//            xmin = x
//            imin = i
//            jmin = j
//         end if
//      end do in                                     // for j = 1 : inner

//      if (flag .eq. one) then
//         x0 = x                                     // save for the next outer iteration
//         call GetH(x0,r) ; r = b-r                  // r = b - matmul(A,x0)
//      else
//         exit out
//      end if

//   end do out                                       // for i = 1 : maxit

//   // returned solution is that with minimum residual
//   if (flag .eq. zero) then
//      relres = normr / n2b
//   else
//      x = xmin
//      iter(1) = imin ; iter(2) = jmin
//      relres = normrmin / n2b
//   end if
//   print*,iter

//   // Cleaning
//   deallocate(V)
//   deallocate(h)
//   deallocate(QT)
//   deallocate(RM)
//   deallocate(f)
//   deallocate(W)
//   deallocate(q)
//   deallocate(y)
//   deallocate(tempv)
}


//_______________________________ DirectSolver ________________________________


DirectSolver :: DirectSolver (boolean store, boolean opt)
    : Solver ()
{
    // Constructor.

    suppress    = opt ;
    storeMatrix = store;

    matrix = NULL;

    if (Mpi_nbProcessors > 1)
	Error("DirectSolver::DirectSolver(opt)",
	      "not implemented in parallel version") ;
}

void DirectSolver :: Solve (SemiDiscreteProblem* prob)
{
    // Solves 'prob', i.e., initializes the unknown field "u" of 'prob' to the
    // solution of the system A u = RHS of 'prob'.

    Field      *unknown, *r ;
    FullMatrix *A ;
    RealVector *sol ;
    int        i, nbDofs ;
    Timer      t_direct("Direct solver") ;

    t_direct.Start() ;

    unknown = prob->GetUnknown() ;
    r       = unknown->FieldGetWork() ;
    nbDofs  = unknown->GetNbDofs() ;

    // 1. right-hand side

    prob->GetRHS(r) ;                     // field form
    sol = new RealVector(nbDofs) ;        // vector form
    r->CopyDofsToVector(sol) ;

    // 2. left-hand side

    if (matrix == NULL){
	printf ("\nGetting the matrix operator in direct solver ... ") ;

	A      = prob->GetMatrixOperator() ;
	matrix = A->Expanded();
	printf ("done\n") ;
    }
    else
	A = matrix->Expanded() ;

    //---------------------------------------------
    // 3. optional b.c. 

    if (suppress) {
	printf("\nsuppress last row and column in DirectSolver\n") ;
	for (i=1 ; i<=nbDofs ; i++) {
	    A->At(i,nbDofs) = 0. ;
	    A->At(nbDofs,i) = 0. ;
	}
	A->At(nbDofs,nbDofs) = 1. ;
	sol->At(nbDofs)      = 0. ;
    }

    // 4. solve

    if (verbosity >= 1) {
	printf("\nOperator LHS: ") ; A->Print() ;
	printf("\nOperator RHS: ") ; sol->Print() ;
    }

    printf("\nEnter Lapack ...") ;
    A->Solve(sol) ;
    printf(" exit Lapack\n") ;

    if (verbosity >= 1) {
	printf("\nsol: ") ; 
	sol->Print() ;
    }

    unknown->CopyDofsFromVector(sol) ;

    delete A ;

    delete sol ;
    unknown->Retrieve(r) ;

    t_direct.Stop() ;
}
