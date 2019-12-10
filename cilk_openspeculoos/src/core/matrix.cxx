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

// matrix.cxx

#include "core/matrix.hxx"
#include "core/options.hxx"
#include "core/timer.hxx"

#include <math.h>
#include <stdio.h>

static real one  = 1. ;
static real zero = 0. ;
static int  ione = 1  ;


//__________________________________ Matrix ___________________________________


void Matrix :: Solve (RealVector*)
{
    // Issues an error message.
 
    ImplementedBySubclasses("Matrix::Solve") ;
}


FullMatrix* Matrix :: Times (FullMatrix* B) 
{
    // Returns a new full matrix, the product of the receiver and 'B'.

# ifdef REQUIRE
    Require("valid dimensions", ncol == B->GetNrow()) ;
    InMethod("DiagonalMatrix::Times(B)") ;
# endif

    FullMatrix* C ;

    C = new FullMatrix(nrow,B->GetNcol()) ;
    C->MatMat(this,B) ;

    return C ;
}


boolean Matrix :: ValidIndices (int i, int j)
{
    // Returns 'true' if the receiver has a coefficient (i,j), else returns
    // 'false'.

    if (i<1 || i>nrow || j<1 || j>ncol) {
	printf ("\nWarning: using index (%d,%d) in Matrix of size (%d,%d)\n",
		i,j,nrow,ncol) ;
	return false ;
    }
    else
	return true ;
}


boolean Matrix :: ValidMatMat (Matrix* A, Matrix* B, char transp)
{
    // Returns True if the dimensions of C (the receiver), 'A' and 'B' make the
    // operation C = A*B (or C = A*B(T), depending on 'transp') legal.

    if (transp == 'N') {
	if (A->nrow == nrow && A->ncol == B->nrow && B->ncol == ncol)
	    return true ;
	else {
	    printf("\nWarning in computing C = A * B : dimensions do not match: C is [%d,%d], A is [%d,%d], B is [%d,%d]\n",nrow,ncol,A->nrow,A->ncol,B->nrow,B->ncol) ;
	    return false ;
	}
    }
    else if (transp == 'T') {
	if (A->nrow == nrow && A->ncol == B->ncol && B->nrow == ncol)
	    return true ;
	else {
	    printf("\nWarning in computing C = A * B(T) : dimensions do not match:C is [%d,%d], A is [%d,%d], B is [%d,%d]\n",nrow,ncol,A->nrow,A->ncol,B->nrow,B->ncol) ;
	    return false ;
	}
    }
    else {
	Error("Matrix::ValidMatMat(A,B,transp)","unexpected value of 'transp'") ;
	return false ;
    }
}


//________________________________ FullMatrix _________________________________


real& FullMatrix :: At (int i, int j)
{
    // Returns a reference to the (i,j)-th coefficient of the receiver.

# ifdef REQUIRE
    Require("valid arguments 'i' and 'j'", ValidIndices(i,j)) ;
    InMethod("FullMatrix::At(i,j)") ;
# endif

    return values[(j-1)*nrow+i-1] ;
}


void FullMatrix :: CopyFrom (FullMatrix* A)
{
    // Initializes the values of the receiver to those of 'A'.

# ifdef REQUIRE
    Require("same dimensions", A->nrow == nrow && A->ncol == ncol) ;
    InMethod("FullMatrix::CopyFrom(A)") ;
# endif

    int n = nrow * ncol ;

    DCOPY(&n, &A->values[0], &ione, &values[0], &ione) ;
}


Matrix* FullMatrix :: Duplicate ()
{
    // Returns a new FullMatrix, the exact copy of the receiver.

    FullMatrix *answer ;

    answer = new FullMatrix(nrow,ncol) ;
    answer->CopyFrom(this) ;
    return answer ;
}


FullMatrix* FullMatrix :: Expanded () 
{
    // Returns a new FullMatrix, the exact copy of the receiver.

    return (FullMatrix*)Duplicate() ;
}

DiagonalMatrix* FullMatrix :: GetEigenpairs (boolean check)
{
    // Computes the eigendecomposition of the receiver A:
    //   A = P lambda inv(P) ,
    // where lambda is the vector of the eigenvalues of A and P is the matrix of
    // the eigenvectors.
    //
    // Overwrites the receiver A with P and returns lambda. The eigenvectors (the
    // columns of A at exit) appear in the same order as the eigenvalues in
    // lambda.
    //
    // Assumes that A is such that its eigenvalues (and thus its eigenvectors)
    // are real; if this is not the case, then:
    //   - if 'check' is 'true' (slower version), an error message is issued,
    //   - if 'check' is 'false' (faster version), lambda and P contain garbage
    //     at exit.
    // If the receiver A is symmetric, this real-eigenvalues condition is met.
    //
    // Warning: this method calls Lapack DGEEV subroutine. DGEEV does not ensure
    //   that the eigenvectors are orthogonal; so do not assume that P(T) is
    //   equal to inv(P).

# ifdef REQUIRE
    Require("square matrix", nrow == ncol) ;
    InMethod("FullMatrix::GetEigenpairs(check)") ;
# endif

    FullMatrix     *vl, *vr ;
    DiagonalMatrix *wr, *wi ;
    RealVector     *work ;
    real           largest ;
    int            i, lwork, info ;
    char           jobvl, jobvr ;
    const real     PREC = 1.e-12 ;

    // 1. Get eigenpairs

    jobvl = 'N' ;                        // do not compute left eigenvectors
    jobvr = 'V' ;                        // compute right eigenvectors
    wr    = new DiagonalMatrix(nrow) ;   // the (real part of the) eigenvalues
    wi    = new DiagonalMatrix(nrow) ;   // the imaginary part of the eigenvalues
    vl    = new FullMatrix(1,nrow) ;     // garbage
    vr    = new FullMatrix(nrow,nrow) ;  // the eigenvectors
    lwork = 4 * nrow ;                   // size of work array
    work  = new RealVector(lwork) ;

    DGEEV (&jobvl, &jobvr, &nrow, values, &nrow, 
	   wr->GetValues(), wi->GetValues(),
	   vl->values, &ione, vr->values, &nrow,
	   work->GetValues(), &lwork, &info) ;

    // 2. Perform some checks

    if (info)
	Error("FullMatrix::GetEigenpairs","error code returned from DGEEV") ;

    if (check) {
	// get largest eigenvalue (real part only); hope not all are 0
	largest = ZERO ;
	for (i=1 ; i<=nrow ; i++)
	    largest = max(largest,fabs(wr->At(i,i))) ;
	// check that all imaginary part of eigenvalues are very small
	for (i=1 ; i<=nrow ; i++)
	    if (fabs(wi->At(i,i)) > PREC*largest)
		Error("FullMatrix::GetEigenpairs","nonreal eigenvalue") ;
    }

    // 3. Overwrite A (which contains garbage now) with its eigenvectors P

    CopyFrom(vr) ;

    delete wi ;
    delete vl ;
    delete vr ;
    delete work ;

    return wr ;
}


void FullMatrix :: Inverse ()
{
    // Overwrites the receiver with its inverse.

# ifdef REQUIRE
    Require("square matrix", nrow == ncol) ;
    InMethod("FullMatrix::Inverse()") ;
# endif

    real *work ;
    int  *ipiv ;
    int  info ;

    if (nrow != 0) {
	ipiv = new int[nrow] ;
	work = new real[nrow] ;
    }
    else {
	ipiv = NULL ;
	work = NULL ;
    }

    // 1. LU factorization

    DGETRF(&nrow, &ncol, &values[0], &nrow, &ipiv[0], &info) ;

# ifdef CHECK
    if (info < 0)
	Error("FullMatrix::Inverse()","invalid value in DGETRF") ;
    else if (info > 0)
	Error("FullMatrix::Inverse()","singular matrix detected in DGETRF") ;
# endif

    // 2. Solves the system to get the inverse

    DGETRI(&nrow, &values[0], &nrow, &ipiv[0], &work[0], &nrow, &info) ;

# ifdef CHECK
    if (info < 0)
	Error("FullMatrix::Inverse()","invalid value in DGETRI") ;
    else if (info > 0)
	Error("FullMatrix::Inverse()","singular matrix detected in DGETRI") ;
# endif

    delete ipiv ;
    delete work ; 
}


void FullMatrix :: MatMat (FullMatrix* A, Matrix* B)
{
    // Initializes the receiver to A*B.
    // Calls DGEMM (BLAS 3).

# ifdef REQUIRE
    Require("Valid matrix dimensions", ValidMatMat(A,B,'N')) ;
    InMethod("FullMatrix::MatMat(FullMatrix* A, Matrix* B)") ;
# endif

    FullMatrix *BB ;
    real       *valuesB ;
    char       trans ;
    int        j, k, size ;

    // E.P.
    switch (B->GetType()) {

	case IDENTITY_MT : 
	    size = nrow * ncol ;
	    //      fprintf(stdout," void FullMatrix :: MatMat (FullMatrix* A, Matrix* B) Dcopy size = %d \n",size);
	    DCOPY(&size, &A->values[0], &ione, &values[0], &ione) ;
	    break ;

	case FULL_MT :
	    trans = 'N' ;
	    BB    = (FullMatrix*)B ;

	    //      fprintf(stdout," void FullMatrix :: MatMat (FullMatrix* A, Matrix* B)FULL_MT DGEMM  size = %d \n",nrow * ncol);

	    DGEMM(&trans, &trans, &A->nrow, &BB->ncol, &A->ncol, &one,
		  &A->values[0], &A->nrow, &BB->values[0], &BB->nrow, &zero, 
		  &values[0], &nrow) ;

	    break ;

	case DIAGONAL_MT :
	    size    = nrow * ncol ;
	    valuesB = B->GetValues() ;

	    //      fprintf(stdout," void FullMatrix :: MatMat (FullMatrix* A, Matrix* B) DIAGONAL_MT Dcopy size = %d \n",size);

	    DCOPY(&size, &A->values[0], &ione, &values[0], &ione) ;
	    k = 0 ;
	    for (j=0 ; j<ncol ; j++) {
		DSCAL(&nrow, &valuesB[j], &values[k], &ione) ;
		k += nrow ;
	    }
	    break ;

	default :
	    Error("FullMatrix::MatMat(FullMatrix*,Matrix*)","Wrong matrix type") ;
    }

}


void FullMatrix :: MatMat (Matrix* A, FullMatrix* B)
{
    // Initializes the receiver to A*B.
    // Calls DGEMM (BLAS 3).

# ifdef REQUIRE
    Require("Valid matrix dimensions", ValidMatMat(A,B,'N')) ;
    InMethod("FullMatrix::MatMat(Matrix* A, FullMatrix* B)") ;
# endif

    FullMatrix *AA ;
    real       *valuesA ;
    char       trans ;
    int        i, size ;
    // E.P.


    //  fprintf(stdout," \n ");

    switch (A->GetType()) {

	case IDENTITY_MT : 
	    size = nrow * ncol ;
	    //      fprintf(stdout," void FullMatrix :: MatMat (Matrix* A, Matrix* B) Dcopy size = %d \n",size);
	    DCOPY(&size, &B->values[0], &ione, &values[0], &ione) ;
	    break ;

	case FULL_MT :
	    trans = 'N' ;
	    AA    = (FullMatrix*)A ;

	    //      fprintf(stderr," void FullMatrix :: MatMat (Matrix* A, Matrix* B)FULL_MT DGEMM n2 = %d \n",AA->ncol );

	    DGEMM(&trans, &trans, &AA->nrow, &B->ncol, &AA->ncol, &one,
		  &AA->values[0], &AA->nrow, &B->values[0], &B->nrow,
		  &zero, &values[0], &nrow) ;
	    break ;

	case DIAGONAL_MT :
	    size    = nrow * ncol ;
	    valuesA = A->GetValues() ;

	    //      fprintf(stdout," void FullMatrix :: MatMat (Matrix* A, Matrix* B) DIAGONAL_MT Dcopy size = %d \n",size);

	    DCOPY(&size, &B->values[0], &ione, &values[0], &ione) ;
	    for (i=0 ; i<nrow ; i++)
		DSCAL(&ncol, &valuesA[i], &values[i], &nrow) ;
	    break ;

	default :
	    Error("FullMatrix::MatMat(Matrix*,FullMatrix*)","Wrong matrix type") ;
    }

}

void FullMatrix :: MatMat (FullMatrix* A, FullMatrix* B)
{
    // Initializes the receiver to A*B.
    // Calls DGEMM (BLAS 3).

# ifdef REQUIRE
    Require("Valid matrix dimensions", ValidMatMat(A,B,'N')) ;
    InMethod("FullMatrix::MatMat(Matrix* A, FullMatrix* B)") ;
# endif

    FullMatrix *AA ;
    real       *valuesA ;
    char       trans ;
    int        i, size ;
    // E.P.


    //  fprintf(stdout," \n ");

    switch (A->GetType()) {

	case IDENTITY_MT : 
	    size = nrow * ncol ;
	    //      fprintf(stdout," void FullMatrix :: MatMat (Matrix* A, Matrix* B) Dcopy size = %d \n",size);
	    DCOPY(&size, &B->values[0], &ione, &values[0], &ione) ;
	    break ;

	case FULL_MT :
	    trans = 'N' ;
	    AA    = (FullMatrix*)A ;

	    //      fprintf(stderr," void FullMatrix :: MatMat (Matrix* A, Matrix* B)FULL_MT DGEMM n2 = %d \n",AA->ncol );

	    DGEMM(&trans, &trans, &AA->nrow, &B->ncol, &AA->ncol, &one,
		  &AA->values[0], &AA->nrow, &B->values[0], &B->nrow,
		  &zero, &values[0], &nrow) ;
	    break ;

	case DIAGONAL_MT :
	    size    = nrow * ncol ;
	    valuesA = A->GetValues() ;

	    //      fprintf(stdout," void FullMatrix :: MatMat (Matrix* A, Matrix* B) DIAGONAL_MT Dcopy size = %d \n",size);

	    DCOPY(&size, &B->values[0], &ione, &values[0], &ione) ;
	    for (i=0 ; i<nrow ; i++)
		DSCAL(&ncol, &valuesA[i], &values[i], &nrow) ;
	    break ;

	default :
	    Error("FullMatrix::MatMat(Matrix*,FullMatrix*)","Wrong matrix type") ;
    }

}


void FullMatrix :: MatTransMat (FullMatrix* A, Matrix* B)
{
    // Initializes the receiver to A * B(T).
    // Calls DGEMM (BLAS 3).

# ifdef REQUIRE
    Require("Valid matrix dimensions", ValidMatMat(A,B,'T')) ;
    InMethod("FullMatrix::MatTransMat(A,B)") ;
# endif

    FullMatrix *BB ;
    real       *valuesB ;
    int        j, k, size ;
    char       transA, transB ;

    // E.P.


    switch (B->GetType()) {

	case IDENTITY_MT : 
	    size = nrow * ncol ;
	    //      fprintf(stdout," void FullMatrix :: MatTransMat (FullMatrix* A, Matrix* B) Dcopy size = %d \n",size);
	    DCOPY(&size, &A->values[0], &ione, &values[0], &ione) ;
	    break ;

	case FULL_MT :
	    transA = 'N' ;
	    transB = 'T' ;
	    BB     = (FullMatrix*)B ;

	    //      fprintf(stdout," void FullMatrix :: MatTransMat (FullMatrix* A, Matrix* B)FULL_MT DGEMM  size = %d \n",nrow * ncol);
	    DGEMM(&transA, &transB, &A->nrow, &BB->nrow, &A->ncol, &one,
		  &A->values[0], &A->nrow, &BB->values[0], &BB->nrow, &zero,
		  &values[0], &nrow) ;
	    break ;

	case DIAGONAL_MT :
	    size    = nrow * ncol ;
	    valuesB = B->GetValues() ;

	    //      fprintf(stdout," void FullMatrix :: MatTransMat (FullMatrix* A, Matrix* B) DIAGONAL_MT Dcopy size = %d \n",size);

	    DCOPY(&size, &A->values[0], &ione, &values[0], &ione) ;
	    k = 0 ;
	    for (j=0 ; j<ncol ; j++) {
		DSCAL(&nrow, &valuesB[j], &values[k], &ione) ;
		k += nrow ;
	    }
	    break ;

	default :
	    Error("FullMatrix::MatTransMat(A,B)","Wrong matrix type") ;
    }

}


void FullMatrix :: MatTransVec (RealVector* x, RealVector* y)
{
    // Initializes y to tA * x, where 'A' stands for the receiver.
    // Calls DGEMV (BLAS 2).

# ifdef REQUIRE
    Require("Valid dimensions", nrow == x->GetSize() && ncol == y->GetSize()) ;
    InMethod("FullMatrix::MatTransVec(x,y)") ;
# endif

    char trans = 'T' ;

    DGEMV(&trans, &nrow, &ncol, &one, &values[0], &nrow,
	  x->GetValues(), &ione, &zero, y->GetValues(), &ione) ;
}


void FullMatrix :: MatVec (RealVector* x, RealVector* y)
{
    // Initializes y to A * x, where 'A' stands for the receiver.
    // Calls DGEMV (BLAS 2).

# ifdef REQUIRE
    Require("Valid dimensions", ncol == x->GetSize() && nrow == y->GetSize()) ;
    InMethod("FullMatrix::Mat(x,y)") ;
# endif

    char trans = 'N' ;
 
    DGEMV(&trans, &nrow, &ncol, &one, &values[0], &nrow,
	  x->GetValues(), &ione, &zero, y->GetValues(), &ione) ;
}


void FullMatrix :: Multiply (real alpha)
{
    // A = alpha A.
    // Multiplies every coefficient of the receiver by 'alpha'.

    int size = nrow*ncol ;
    int ione = IONE ;

    // test due to a bug generated by PURIFY or BLAS!
# ifndef PURIFY
    DSCAL(&size, &alpha, &values[0], &ione) ;
# else
    for (int i=0 ; i<size ; i++)
	values[i] *= alpha ;
# endif
}


void FullMatrix :: Print ()
{
    // Prints the receiver on standard output.

    int i, j, max_row, max_col ;

    max_row = 20 ;
    max_col = 20 ;

    printf("FullMatrix of size %d by %d:\n",nrow,ncol) ;

    if (ncol > max_col)
	printf ("   only the first %d columns are printed\n",max_col) ;

    if (nrow > max_row)
	printf ("   only the first %d rows are printed\n",max_row) ;

    for (i=1 ; i<=nrow; ++i) {
	printf (" ") ;
	for (j=1 ; j<=ncol && j<=max_col; ++j)
//      printf ("% .2e ",At(i,j)) ;
	    printf ("%le ",At(i,j)) ;
	printf ("\n") ;
    }
}

void FullMatrix :: PrintMatlab() 
{

    int i, j ; 
    FILE *ii_f, *jj_f, *aa_f ;

    ii_f = fopen("ii.mat", "w");
    jj_f = fopen("jj.mat", "w");
    aa_f = fopen("aa.mat", "w");
  
    if( ii_f == NULL || jj_f == NULL || aa_f == NULL ) {
	printf(" Not enougth memory for aa ii jj \n") ; 
	exit(-1) ; 
    }
  
    for (i=1 ; i<=nrow ; i++) {
	for (j=1 ; j<=ncol ; j++) {
	    fprintf (ii_f, "%d\n", i);
	    fprintf (jj_f, "%d\n", j);
	    fprintf (aa_f, "%le\n", At(i,j));
	}
    }
  
    fclose(ii_f);
    fclose(jj_f);
    fclose(aa_f);
  
}

void FullMatrix :: PrintColumn (int j)
{
    // Prints on standard output the 'j'-th column of the receiver.

# ifdef REQUIRE
    Require("valid argument 'j'", j>=1 || j<=ncol) ;
    InMethod("FullMatrix::PrintColumn(j)") ;
# endif

    int i ;
    const int PER_LINE = 6 ;   // number of values printed on every line

    printf("Column %d of FullMatrix of size %d by %d:\n",j,nrow,ncol) ;

    for (i=1  ; i<=nrow ; i++) {
	printf (" % .4e",At(i,j)) ;
	if (i%PER_LINE==0 && i!=nrow)
	    printf ("\n  ") ;
    }

    printf ("\n") ;
}


void FullMatrix :: Solve (RealVector* x)
{
    // Overwrites 'x' with 'inv(A) * x', where 'A' stands for the receiver.

# ifdef REQUIRE
    Require("square matrix", nrow == ncol) ;
    Require("valid dimensions", x->GetSize() == ncol) ;
    InMethod("FullMatrix::Solve(x)") ;
# endif

    int *ipiv, info ;

    if (nrow)
	ipiv = new int[nrow] ;
    else
	ipiv = NULL ;

    DGESV(&ncol, &ione, &values[0], &nrow, &ipiv[0], 
	  x->GetValues(), &nrow, &info);

# ifdef CHECK
    if (info < 0)
	Error("FullMatrix::Solve(x)","illegal value of some coefficient") ;
    if (info > 0)
	Error("FullMatrix::Solve(x)","singular matrix") ;
# endif   

    delete ipiv ;
}


void FullMatrix :: Solve (FullMatrix* X, FullMatrix* B)
{
    // Initializes 'X' to 'inv(A) * B', where 'A' stands for the receiver.

# ifdef REQUIRE
    Require("square matrix", nrow == ncol) ;
    Require("A and X are consistent", X->nrow = ncol) ;
    Require("X and B are consistent", X->nrow == B->nrow && X->ncol == B->ncol) ;
    InMethod("FullMatrix::Solve(X,B)") ;
# endif

    int *ipiv, info ;

    if (nrow)
	ipiv = new int[nrow] ;
    else
	ipiv = NULL ;

    X->CopyFrom(B) ;

    DGESV(&ncol, &X->ncol, &values[0], &nrow, &ipiv[0], 
	  &X->values[0], &nrow, &info);

# ifdef CHECK
    if (info < 0)
	Error("FullMatrix::Solve(X,B)","illegal value of some coefficient") ;
    if (info > 0)
	Error("FullMatrix::Solve(X,B)","singular matrix") ;
# endif   

    delete ipiv ;
}


boolean FullMatrix :: TestSymmetry (boolean print)
{
    // Returns 'true' if the receiver is symmetric, else returns 'false'.
    // If 'debug' is set to 'true', prints information on the location of the
    // non-symmetric coefficients, if any.

    int     i, j ;
    boolean answer ;

    if (nrow != ncol)
	return false ;

    answer = true ;
    for (j=1 ; j<=ncol ; j++)
	for (i=j+1 ; i<=nrow ; i++) 
	    if (fabs(At(i,j)-At(j,i)) > 1.e-7) {
		answer = false ;
		if (print) 
		    printf("not symmetric in %i %i : %e %e\n",i,j,At(i,j),At(j,i)) ;
	    }

    return answer ;
}


void FullMatrix :: TransMatMat (Matrix* A, FullMatrix* B)
{
    // Initializes the receiver to tA * B.
    // Calls DGEMM (BLAS 3).

# ifdef REQUIRE
    Require("Valid matrix dimensions",
	    A->GetNcol() == nrow && A->GetNrow() == B->nrow && B->ncol == ncol);
    InMethod("FullMatrix::TransMatMat(A,B)") ;
# endif

    FullMatrix *AA ;
    char       transA, transB ;
    int        size ;
    // E.P.



    switch (A->GetType()) {
	case IDENTITY_MT : 
	    size = nrow * ncol ;
	    DCOPY(&size, &B->values[0], &ione, &values[0], &ione) ;
	    break ;
	case FULL_MT :
	    transA = 'T' ;
	    transB = 'N' ;
	    AA     = (FullMatrix*)A ;
	    DGEMM(&transA, &transB, &AA->ncol, &B->ncol, &AA->nrow, &one,
		  &AA->values[0], &AA->nrow, &B->values[0], &B->nrow, &zero, 
		  &values[0], &nrow) ;
	    break ;
	default :
	    Error("FullMatrix::TransMatMat","Wrong matrix type") ;
    }

}


//______________________________ FullMatrixLU _________________________________


FullMatrixLU :: FullMatrixLU (int n, int m)
    : FullMatrix(n,m)
{
    // Constructor.

# ifdef REQUIRE
    Require("square matrix", nrow == ncol) ;
    InMethod("FullMatrixLU::FullMatrixLU(n,m)") ;
# endif

    isFactorized = false ;
    ipiv         = NULL ;
}


FullMatrixLU :: ~FullMatrixLU ()
{
    // Destructor.

    delete ipiv ;
}


void FullMatrixLU :: Factorize ()
{
    // Overwrites the receiver with its LU decomposition.

# ifdef REQUIRE
    Require("matrix is not factorized", isFactorized == false) ;
    InMethod("FullMatrixLU::Factorize(x)") ;
# endif
  
    int info ;

    if (nrow)
	ipiv = new int[nrow] ;

    DGETRF(&nrow, &ncol, &values[0], &nrow, &ipiv[0], &info) ;

# ifdef CHECK
    if (info < 0)
	Error("FullMatrixLU::Factorize()","invalid value in DGETRF") ;
    else if (info > 0)
	Error("FullMatrixLU::Factorize()","singular matrix detected in DGETRF") ;
# endif

    isFactorized = true ;
}


void FullMatrixLU :: Solve (RealVector* x)
{
    // Overwrites 'x' with inv(A) * x, where 'A' stands for the receiver.

# ifdef REQUIRE
    Require("valid dimensions", x->GetSize() == ncol) ;
    InMethod("FullMatrixLU::Solve(x)") ;
# endif

    int  info ;
    char trans ;

    if (! isFactorized)
	Factorize() ;

    trans = 'N' ;

    DGETRS(&trans, &ncol, &ione, &values[0], &nrow, &ipiv[0], 
	   x->GetValues(), &nrow, &info);

# ifdef CHECK
    if (info < 0)
	Error("FullMatrixLU::Solve(x)","illegal value of some coefficient") ;
# endif   
}


//_______________________________ DiagonalMatrix ______________________________


DiagonalMatrix :: DiagonalMatrix (int n)
    : Matrix (n,n)
{
    // Constructor. Initializes the receiver as a diagonal matrix of size (n,n).

    if (n)
	values = new real[n] ;
    else
	values = NULL ;
}


DiagonalMatrix :: ~DiagonalMatrix ()
{
    // Destructor.

      delete values ;
}


DiagonalMatrix :: DiagonalMatrix (RealVector* x)
    : Matrix (x->GetSize(),x->GetSize())
{
    // Constructor. Initializes the receiver to a diagonal matrix containing the
    // same values as 'x'.
 
    int i ;

    if (nrow)
	values = new real[nrow] ;
    else
	values = NULL ;

    for (i=1 ; i<=nrow ; i++)
	At(i,i) = x->At(i) ;
}


real& DiagonalMatrix :: At (int i, int j)
{
    // Returns the (i,j)-th coefficient of the receiver.

# ifdef REQUIRE
    Require("valid arguments (i,j)", i >= 1 && i <= nrow && j == i) ;
    InMethod("DiagonalMatrix::At(i,j)") ;
# endif

    int foo = j ;        // just to avoid a compilation warning!

    return values[i-1] ;
}


Matrix* DiagonalMatrix :: Duplicate ()
{
    // Returns a new DiagonalMatrix, the exact copy of the receiver.

    DiagonalMatrix *answer ;
    int            i ;

    answer = new DiagonalMatrix(nrow) ;
    for (i=1 ; i<=nrow ; i++)
	answer->At(i,i) = At(i,i) ;

    return answer ;
}


FullMatrix* DiagonalMatrix :: Expanded () 
{
    // Returns a new FullMatrix, the receiver in full form.

    FullMatrix *answer ;
    int        i ;

    answer = new FullMatrix(nrow,nrow) ;
    for (i=1 ; i<=nrow ; i++)
	answer->At(i,i) = At(i,i) ;

    return answer ;
}


real DiagonalMatrix :: GetInfiniteNorm () 
{
    // Returns the infinite norm (= max norm) of the receiver.

    real answer ;
    int  i ;

    answer = ZERO ;
    for (i=0 ; i<nrow ; i++)
	answer = max(answer,fabs(values[i])) ;

    return answer ;
}


boolean DiagonalMatrix :: HasZeroes ()
{
    // Returns 'true' if one of coefficients of the receiver is zero, else 
    // returns 'false'.

    int i ;

    for (i=0 ; i<nrow ; i++)
	if (values[i] == ZERO)
	    return true ;

    return false ;
}


void DiagonalMatrix :: Inverse ()
{
    // Replaces every coefficient of the receiver by its inverse.

# ifdef REQUIRE
    Require("has no zeroes", ! HasZeroes()) ;
    InMethod("DiagonalMatrix::Inverse()") ;
# endif

    int i ;

    for (i=0 ; i<nrow ; i++)
	values[i] = ONE/values[i] ;
}


void DiagonalMatrix :: MatVec (RealVector* x, RealVector* y)
{
    // Initializes y to A * x, where 'A' stands for the receiver.

# ifdef REQUIRE
    Require("Valid dimensions", ncol == x->GetSize() && nrow == y->GetSize()) ;
    InMethod("DiagonalMatrix::Mat(x,y)") ;
# endif

    real *valX, *valY ;
    int  i ;

    valX = x->GetValues() ;
    valY = y->GetValues() ;

    for (i=0 ; i<nrow ; i++)
	valY[i] = values[i] * valX[i] ;
} 


void DiagonalMatrix :: MatTransVec (RealVector*, RealVector*)
{
    // Not implemented yet.

    NotImplemented("DiagonalMatrix::MatTransVec") ;
}


void DiagonalMatrix :: Print ()
{
    // Prints the receiver on standard output.

    int i ;
    const int PER_LINE = 6 ;        // number of values printed on every line

    printf("DiagonalMatrix of size %d:\n",nrow);

    for (i=1  ; i<=nrow ; i++) {
//    printf (" % .4e",At(i,i)) ;
	printf (" %le",At(i,i)) ;
	if (i%PER_LINE==0 && i!=nrow)
	    printf ("\n") ;
    }

    printf ("\n") ;
}


//________________________________ IdentityMatrix _____________________________


IdentityMatrix :: IdentityMatrix (int n)
    : DiagonalMatrix (n)
{
    // Constructor. Initializes the receiver to the identity matrix of rank 'n'.

    int i ;

    for (i=0 ; i<nrow ; i++)
	values[i] = ONE ;
}


//________________________________ TensorMatrix _______________________________

Pool* TensorMatrix :: pool = new Pool("TensorMatrix",sizeof(TensorMatrix),1) ;
// Initializes the class attribute 'pool'.


FullMatrix* TensorMatrix :: auxiliaryMatrix1 = new FullMatrix(20000,1) ;
// Initializes the static attribute 'auxiliaryMatrix1'.


FullMatrix* TensorMatrix :: auxiliaryMatrix2 = new FullMatrix(20000,1) ;
// Initializes the static attribute 'auxiliaryMatrix1'.


void TensorMatrix :: DeleteMatrices ()
{
    // Deletes the directional matrices of the receiver. By default, when a 
    // tensor matrix is deleted, its directional matrices are not deleted.

    delete matrix1 ;
    delete matrix2 ;
    delete matrix3 ;
}


FullMatrix* TensorMatrix :: Expanded ()
{
    // A debugging tool.
    // Returns a new matrix, the full matrix equivalent of the receiver.

    FullMatrix *answer, *A1, *A2, *A3 ;
    int        i1, i2, i3, j1, j2, j3, n1, n2, n3, m1, m2, m3, ii, jj ;

    answer = new FullMatrix(GetNrow(),GetNcol()) ;

    if (dimension == 1) {
	A1 = matrix1->Expanded() ;    // because message matrix1->At(1,2) is not
	n1 = A1->GetNrow() ;          // valid if matrix1 is not a FullMatrix
	m1 = A1->GetNcol() ;
	for (i1=1; i1 <= n1 ; i1++)
	    for (j1=1 ; j1 <= m1 ; j1++)
		answer->At(i1,j1) = A1->At(i1,j1) ;
	delete A1 ;
    }
    else if (dimension == 2) {
	A1 = matrix1->Expanded() ;
	A2 = matrix2->Expanded() ;
	n1 = A1->GetNrow() ;
	n2 = A2->GetNrow() ;
	m1 = A1->GetNcol() ;
	m2 = A2->GetNcol() ;
	for (i1=1; i1 <= n1 ; i1++)
	    for (i2=1 ; i2 <= n2 ; i2++)
		for (j1=1 ; j1 <= m1 ; j1++)
		    for (j2=1 ; j2 <= m2 ; j2++) {
			ii = (i2-1)*n1 + i1 ;
			jj = (j2-1)*m1 + j1 ;
			answer->At(ii,jj) = A1->At(i1,j1) * A2->At(i2,j2) ;
		    }
	delete A1 ;
	delete A2 ;
    }

    else {
	A1 = matrix1->Expanded() ;
	A2 = matrix2->Expanded() ;
	A3 = matrix3->Expanded() ;
	n1 = A1->GetNrow() ;
	n2 = A2->GetNrow() ;
	n3 = A3->GetNrow() ;
	m1 = A1->GetNcol() ;
	m2 = A2->GetNcol() ;
	m3 = A3->GetNcol() ;
	for (i1=1; i1 <= n1 ; i1++)
	    for (i2=1 ; i2 <= n2 ; i2++)
		for (i3=1 ; i3 <= n3 ; i3++)
		    for (j1=1 ; j1 <= m1 ; j1++)
			for (j2=1 ; j2 <= m2 ; j2++)
			    for (j3=1 ; j3 <= m3 ; j3++) {
				ii = (i3-1)*n1*n2 + (i2-1)*n1 + i1 ;
				jj = (j3-1)*m1*m2 + (j2-1)*m1 + j1 ;
				answer->At(ii,jj) = A1->At(i1,j1) * A2->At(i2,j2) 
				    * A3->At(i3,j3) ;
			    }
	delete A1 ;
	delete A2 ;
	delete A3 ;
    }

    return answer ;
}


Matrix* TensorMatrix :: GetMatrix (int i)
{
    // Returns the 'i'-th directional matrix of the receiver.

# ifdef REQUIRE
    Require("valid argument 'i'", i >= 1 && i <= dimension) ;
    InMethod("TensorMatrix::GetMatrix(i)") ;
# endif

    if (i == 1)
	return matrix1 ;
    else if (i == 2)
	return matrix2 ;
    else
	return matrix3 ;
}


int TensorMatrix :: GetNcol ()
{
    // Returns the number of columns of the matrix represented by the receiver.

    if (dimension == 1)
	return matrix1->GetNcol() ;

    else if (dimension == 2)
	return matrix1->GetNcol() * matrix2->GetNcol() ;

    else
	return matrix1->GetNcol() * matrix2->GetNcol() * matrix3->GetNcol() ;
}


int TensorMatrix :: GetNrow ()
{
    // Returns the number of rows of the matrix represented by the receiver.

    if (dimension == 1)
	return matrix1->GetNrow() ;

    else if (dimension == 2)
	return matrix1->GetNrow() * matrix2->GetNrow() ;

    else
	return matrix1->GetNrow() * matrix2->GetNrow() * matrix3->GetNrow() ;
}


boolean TensorMatrix :: IsConsistent ()
{
    // Returns True if the number of directional matrices of the receiver is
    // equal to the dimension of the receiver, else returns False.

    if (dimension == 1)
	return (matrix1 && !matrix2 && !matrix3) ;

    else if (dimension == 2)
	return (matrix1 && matrix2 && !matrix3) ;

    else
	return (matrix1 && matrix2 && matrix3) ;
}


DiagonalMatrix* TensorMatrix :: Inverse (boolean modifyZeroEigenvalues)
{
    // Computes the explicit inversion of the receiver, by means of Lynch's fast
    // diagonalization technique.
    //
    // Assumes, exceptionally, that the receiver stores the following tensor
    // matrix:
    //   T = I o I o A1 + I o A2 o I + A3 o I o I  (in 3D),
    // instead of T = A3 o A2 o A1, where "o" designates the tensor product.
    //
    // The explicit inverse of T is given by
    //   inv(T) = (P3 o P2 o P1) inv(lambda) (inv(P3) o inv(P2) o inv(P1)),
    // where:
    //   . lambda = I o I o lambda1 + I o lambda2 o I + lambda3 o I o I,
    //   . inv(P1) A1 P1 = lambda1, (and likewise for P2 and P3),
    //   . lambda1 is the diagonal matrix of the eigenvalues of A1,
    //   . P1 is the matrix of the eigenvectors of A1; inv(P1) = P1(transp).
    //
    // Inverting lambda fails if each directional eigenvalue matrix lambda_i has
    // a null eigenvalue, because it implies that lambda has one null eigenvalue.
    // In that case:
    //   . if 'modifyZeroEigenvalues' is 'true', the (infinite) inverse of the
    //     null eigenvalue is replaced by zero, which means that the contribution
    //     of the parasitic mode is omitted;
    //   . if 'modifyZeroEigenvalues' is 'false', an error message is issued.
    //
    // At method exit,
    //   . P1,P2,P3 overwrite A1,A2,A3, respectively;
    //   . the answer is a diagonal matrix of size np1*np2*np3 (in 3D), which
    //     contains the inverse of lambda. Note that inv(lambda), unlike lambda,
    //     has no tensor form. 
    //
    // Reference: R.E. Lynch, J.R. Rice, D.H. Thomas, Direct Simulation of 
    //            Partial Difference Equations by Tensor Products Methods,
    //            Numerische Mathematik, Vol. 6, pp 185-199, 1964.
    // See also:  W. Couzy and M.O. Deville, Spectral-Element Preconditioners for
    //            the Uzawa Pressure Operator Applied to Incompressible Flows,
    //            J. of Scientific Computing, Vol. 9, No. 2, pp 107-122, 1994.
    //            (also in Couzy's thesis report, p 96)

    TensorMatrix   *lambda ;
    FullMatrix     *a ;
    DiagonalMatrix *lambda_i, *answer ;
    RealVector     *one, *temp ;
    real           largest, eigenval ;
    int            i, j, n ;
    boolean        check ;
    const real     SMALL=1.e-10 ;

    // 1. Compute the eigenpairs (lambda_i, P_i) of A_i
    //    lambda_i is a diagonal matrix, P_i overwrites A_i

    lambda = new TensorMatrix(dimension) ;

    for (i=1 ; i<=dimension ; i++) {

	a = (FullMatrix*)GetMatrix(i) ;

#   ifdef CHECK
	check = true ;
	if (a->GetType() != FULL_MT)
	    Error("TensorMatrix::Inverse()","A(i) should be a full matrix") ;
#   else
	check = false ;
#   endif

	lambda_i = a->GetEigenpairs(check) ;        // P_i ovewrites A_i
	lambda->SetMatrix(i,lambda_i) ;
    }

    // 2. Get the explicit form of lambda, i.e., expand
    //    lambda = I o I o lambda1  +  I o lambda2 o I  +  lambda3 o I o I

    n    = lambda->GetNrow() ;
    one  = new RealVector(n) ;
    temp = new RealVector(n) ;

    one->SetValues(ONE) ;
    lambda->MatVec(one,temp) ;
    answer = new DiagonalMatrix(temp) ;
    
    // 3. Invert lambda

    if (modifyZeroEigenvalues) {
	// detect the null eigenvalues of lambda (not of the lambda_i's) and
	// replace their inverse by 0
	largest = answer->GetInfiniteNorm() ;
	for (j=1 ; j<=n ; j++) {
	    eigenval = answer->At(j,j) ;
	    if (fabs(eigenval) < SMALL * largest)        // null eigenvalue
		answer->At(j,j) = ZERO ;
	    else                                         // non-null eigenvalue
		answer->At(j,j) = ONE / eigenval ; 
	}
    }
    else
	answer->Inverse() ;

    // 4. Cleaning

    lambda->DeleteMatrices() ;
    delete lambda ;
    delete one ;
    delete temp ;
    
    return answer ;
}


void TensorMatrix :: MatTransVec (RealVector* x, RealVector* y)
{
    // Initializes 'y' to A * x, where 'A' stands for the receiver.
    //
    // For example : in 3D, Y = (A3(T) o A2(T) o A1(T)) * x , where: 
    //   - "o" denotes the tensor product,
    //   - (T) denotes a matrix transposition,
    //   - the dimensions are :  Y (I*J*K), X (L*M*N),
    //                           A1 (L*I), A2 (M*J), A3 (N*K).

# ifdef REQUIRE
    Require("x has correct dimension", x->GetSize() == GetNrow()) ;
    Require("y has correct dimension", y->GetSize() == GetNcol()) ;
    Require("valid dimension", dimension >= 1 && dimension <= 3) ;
    Require("is consistent", IsConsistent()) ;
    InMethod("TensorMatrix::MatTransVec(x,y)") ;
# endif

    Matrix     *A1, *A2, *A3 ;
    FullMatrix *work1, *work2, *work3, *Xmat ;
    int        k, size, work2_nrow ;
    int        I=1, J=1, K=1, L=1, M=1, N=1 ;
    char       transA, transB ;

    // E.P.



    transA = 'N' ;
    transB = 'N' ;

    A1  = matrix1 ;
    I   = A1->GetNcol() ;
    L   = A1->GetNrow() ;

    if (dimension > 1) {
	A2 = matrix2 ;
	J  = A2->GetNcol() ;
	M  = A2->GetNrow() ;
	if (dimension > 2) {
	    A3 = matrix3 ;
	    K  = A3->GetNcol() ;
	    N  = A3->GetNrow() ;
	}
    }

    // 1. A1 (first dimension)

    if (dimension == 1)
	A1->MatTransVec(x,y) ;

    else {
	// get work matrices
	work1 = auxiliaryMatrix1 ;
	work1->SetDimensions(I,M*N) ;
	Xmat  = new FullMatrix(x,L,M*N) ;            // aliased
	work1->TransMatMat(A1,Xmat) ;

	// 2. A2 (second dimension)

	if (dimension == 2)
	    // create an 'alias' matrix to vector y
	    work2 = new FullMatrix(y,I,J) ;
	else {
	    // get a temporary storage
	    work2 = auxiliaryMatrix2 ;
	    work2->SetDimensions(I*J,N) ;
	}
	work2_nrow = I ;

	work1->SetDimensions(I*M,N) ;
	switch (A2->GetType()) {
	    case IDENTITY_MT : 
		size = I*M*N ; 
		DCOPY(&size, work1->GetValues(), &ione, work2->GetValues(), &ione) ;
		break ;
	    case FULL_MT :
		for (k=1; k<=N; k++)
		    DGEMM(&transA, &transB, &I, &J, &M, &one,
			  work1->GetColumn(k), &I, A2->GetValues(), &M, &zero,
			  work2->GetColumn(k), &work2_nrow) ;
		break ;
	    default :
		Error("TensorMatrix::MatTransVec(x,y)",
		      "Wrong directional matrix type") ;
	}

	// 3. A3 ( third dimension )

	if (dimension == 3) {
	    A3    = matrix3 ;
	    work3 = new FullMatrix(y,I*J,K) ;     // aliased
	    work3->MatMat(work2,A3) ;
	    delete work3 ;
	}

	// 4. Cleaning

	if (dimension == 2)
	    delete work2 ;             // was aliased

	delete Xmat ;                // was aliased
    }
}


void TensorMatrix :: MatVec (RealVector* x, RealVector* y)
{
    // Initializes 'y' to A * x, where 'A' stands for the receiver.
    //
    // For example : in 3D, y = (A3 o A2 o A1) * x , where:
    //   - "o" denotes the tensor product,
    //   - the dimensions are :  Y (L*M*N), X (I*J*K),
    //                           A1 (L*I), A2 (M*J), A3 (N*K).
    //   - the coefficients of X and Y are stored as in elementary field, i.e.,
    //     layer by layer; on every layer, they are stored row by row.

# ifdef REQUIRE
    Require("x has correct dimension", x->GetSize() == GetNcol()) ;
    Require("y has correct dimension", y->GetSize() == GetNrow()) ;
    Require("valid dimension", dimension >= 1 && dimension <= 3) ;
    Require("is consistent", IsConsistent()) ;
    InMethod("TensorMatrix::MatVec(x,y)") ;
# endif

    Matrix     *A1, *A2, *A3 ;
    FullMatrix *work1, *work2, *work3, *Xmat ;
    int        j, k, size, work2_nrow ;
    int        I=1, J=1, K=1, L=1, M=1, N=1 ;
    char       transA, transB ;

    transA = 'N' ;
    transB = 'T' ;

    A1  = matrix1 ;
    I   = A1->GetNcol() ;
    L   = A1->GetNrow() ;

    if (dimension > 1) {
	A2 = matrix2 ;
	J  = A2->GetNcol() ;
	M  = A2->GetNrow() ;
	if (dimension > 2) {
	    A3 = matrix3 ;
	    K  = A3->GetNcol() ;
	    N  = A3->GetNrow() ;
	}
    }
    // 1. A1 (first dimension)

    if (dimension == 1)
	// create an 'alias' matrix to the Y vector
	work1 = new FullMatrix(y,L,1) ;
    else {
	// get a matrix from a free list
	work1 = auxiliaryMatrix1 ;
	work1->SetDimensions(L,J*K) ;
    }

    Xmat = new FullMatrix(x,I,J*K) ;    // aliased
    work1->MatMat(A1,Xmat) ;

    // 2. A2 (second dimension)

    if (dimension > 1) {
	if (dimension == 2)
	    // create an 'alias' matrix to the Y vector
	    work2 = new FullMatrix(y,L,M) ;
	else {
	    // create a storage from the list
	    work2 = auxiliaryMatrix2 ;
	    work2->SetDimensions(L*M,K) ;
	}
	work2_nrow = L ;

	work1->SetDimensions(L*J,K) ;
	switch (A2->GetType()) {
	    case IDENTITY_MT : 
		size = L*J*K ; 
		DCOPY(&size, work1->GetValues(), &ione, work2->GetValues(), &ione) ;
		break ;
	    case FULL_MT :
		for (k=1 ; k<=K ; k++)
		    DGEMM(&transA, &transB, &L, &M, &J, &one,
			  work1->GetColumn(k), &L, A2->GetValues(), &M, &zero,
			  work2->GetColumn(k), &work2_nrow) ;
		break ;
	    case DIAGONAL_MT :
		size = L*J*K ; 
		DCOPY(&size, work1->GetValues(), &ione, work2->GetValues(), &ione) ;
		for (k=1 ; k<=K ; k++)
		    for (j=1 ; j<=J ; j++)
			DSCAL(&L, &A2->At(j,j), work2->GetColumn(k), &ione) ;
		break ;
	    default :
		Error("TensorMatrix::MatVec(x,y)","Wrong directional matrix type") ;
	}

	// 3. A3 (third dimension)

	if (dimension == 3) {
	    work3 = new FullMatrix(y,L*M,N) ;      // aliased
	    work3->MatTransMat(work2,A3) ;
	    delete work3 ;
	}

	// 4. Cleaning

	if (dimension == 2)
	    delete work2 ;           // was aliased
    }

    if (dimension == 1)
	delete work1 ;             // was aliased

    delete (Xmat) ;              // was aliased
}


void TensorMatrix :: Print ()
{
    // Prints the receiver on standard output.

    printf("Tensor matrix of order %d\n",dimension) ;

    printf("matrix 1: ") ;
    if (matrix1)
	matrix1->Print() ;
    else
	printf("NULL\n") ;

    if (dimension >= 2) {
	printf("matrix 2: ") ;
	if (matrix2)
	    matrix2->Print() ;
	else
	    printf("NULL\n") ;
    }

    if (dimension == 3) {
	printf("matrix 3: ") ;
	if (matrix3)
	    matrix3->Print() ;
	else
	    printf("NULL\n") ;
    }
}


void TensorMatrix :: SetMatrix (int i, Matrix* mat)
{
    // Sets the 'i'-th directional matrix of the receiver to 'mat' (no copy).

# ifdef REQUIRE
    Require("valid argument 'i'", i >= 1 && i <= dimension) ;
    InMethod("TensorMatrix::SetMatrix(i,mat)") ;
# endif

    if (i == 1)
	matrix1 = mat ;

    else if (i == 2)
	matrix2 = mat ;
   
    else
	matrix3 = mat ;
}


