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

#ifndef MATRIX
#define MATRIX

/*!
 * \file matrix.hxx
 * \brief Manages matrices
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/vector.hxx"
#include "core/pool.hxx"
#include "core/options.hxx"

class FullMatrix ;
class DiagonalMatrix ;
enum MatrixType {FULL_MT, DIAGONAL_MT, IDENTITY_MT} ;


/*! \class Matrix
   * \brief Base class of matrix classes (except TensorMatrix).
   * 
   * This abstract class provides the common features to matrix classes.
   * implementations. 
   */ 

class Matrix
{

protected :

    int   nrow ;
    int   ncol ;
    real* values ;

    boolean ValidIndices (int i, int j) ;
    boolean ValidMatMat (Matrix* A, Matrix* B, char transp) ;

public :

    // = CONSTRUCTORS

    inline          Matrix (int n, int m) ;
    inline virtual ~Matrix () ;

    // = INQUIRY

    virtual real&       At (int i, int j) = 0 ;
    inline int          GetNcol () ;
    inline int          GetNrow () ;
    virtual MatrixType  GetType () = 0 ;
    inline real*        GetValues () ;
    virtual void        Print () = 0 ;

    // = OPERATIONS

    virtual Matrix*     Duplicate () = 0 ;
    virtual void        Inverse () = 0 ;
    virtual void        MatVec (RealVector* x, RealVector* y) = 0 ;
    virtual void        MatTransVec (RealVector* x, RealVector* y) = 0 ;
    virtual void        Solve (RealVector* x) ;
    FullMatrix*         Times (FullMatrix* B) ;

    // = DEBUGGING

    virtual FullMatrix* Expanded () = 0 ;
} ;


/*! \class FullMatrix
   * \brief Implements rectangular matrices.
   * 
   * The coefficients are stored in <{values}>, columnwise (as in Fortran).
   * 
   * The dimension (nrow,ncol) of a matrix can be modified, provided that
   * the dimension does not exceed <{allocatedSize}>, the size of
   * <{values}>.
   * 
   * A matrix A can be aliased to another matrix B. This implies that
   * attribute <{values}> of A is set to that of B. This allows the user
   * to consider B as A with modified number of rows and columns; for 
   * example, A is 4x4 and B is 8x2. Accordingly, if A is deleted,
   * <{values}> cannot be deleted. The additional attribute <{isAlias}> 
   * holds this information; its default value is 'false'.
   */ 

class FullMatrix : public Matrix
{

protected :

    int     allocatedSize ;
    boolean isAlias ;

public :

    // = CONSTRUCTORS

    inline FullMatrix () ;
    inline FullMatrix (int n,  int m) ;
    inline FullMatrix (RealVector* v, int n, int m) ;
    inline virtual ~FullMatrix () ;

    // = ACCESS TO THE VALUES
 
    real&              At (int i, int j) ;
    inline real*       GetColumn (int j) ;
 
    // = DEFINITIONS

    inline void        DefineAlias (RealVector* x, int n, int m) ;
    inline void        SetDimensions (int n, int m) ;

    // = INQUIRY
 
    inline int         GetAllocatedSize () ;
    inline MatrixType  GetType () ;

    // = MATRIX-MATRIX PRODUCTS

    void               MatMat (Matrix* A, FullMatrix* B) ;
    void               MatMat (FullMatrix* A, Matrix* B) ;
    void               MatMat (FullMatrix* A, FullMatrix* B) ;
    void               MatTransMat (FullMatrix* A, Matrix* B) ;
    void               TransMatMat (Matrix* A, FullMatrix* B) ;

    // = MATRIX-VECTOR PRODUCTS

    void               MatVec (RealVector* x, RealVector* y) ;
    void               MatTransVec (RealVector* x, RealVector* y) ;

    // = OTHER OPERATIONS

    void               CopyFrom (FullMatrix* A) ;
    Matrix*            Duplicate () ;
    DiagonalMatrix*    GetEigenpairs (boolean check=false) ;
    void               Inverse () ;
    void               Multiply (real alpha) ;
    virtual void       Solve (RealVector* x) ;
    void               Solve (FullMatrix* X, FullMatrix* B) ;

    // = DEBUGGING
  
    FullMatrix*        Expanded () ;
    void               Print () ;
    void               PrintColumn (int j) ;
    void               PrintMatlab () ; 
    boolean            TestSymmetry (boolean print=false) ;
   
} ;


/*! \class FullMatrixLU
   * \brief FIXME: no doc
   * 
   * FIXME : [Vincent Keller] No documentation !
   */ 

class FullMatrixLU : public FullMatrix
{

protected :

    boolean isFactorized ;
    int*    ipiv ;

    void Factorize () ;

public :

    // = CONSTRUCTORS

    FullMatrixLU (int n,  int m) ;
    ~FullMatrixLU () ;

    // = SOLVING

    void Solve (RealVector* x) ;
} ;


/*! \class DiagonalMatrix
   * \brief Manages diagonal matrices
   * 
   * Diagonal matrices are stored as arrays of size nrow=ncol.
   */ 
class DiagonalMatrix : public Matrix
{

public :

    // = CONSTRUCTORS

    DiagonalMatrix (int n) ;
    DiagonalMatrix (RealVector* x) ;
    ~DiagonalMatrix () ;

    // OPERATORS

    real&             At (int i, int j) ;

    // INQUIRY

    inline MatrixType GetType () ;
    boolean           HasZeroes () ;
    void              Print () ;

    // = OPERATIONS
 
    Matrix*           Duplicate () ;
    void              MatVec (RealVector* x, RealVector* y) ; 
    void              MatTransVec (RealVector* x, RealVector* y) ;
    real              GetInfiniteNorm () ;
    void              Inverse () ;

    // = DEBUGGING

    FullMatrix*       Expanded () ;
} ;

 
/*! \class IdentityMatrix
   * \brief Manages identity matrices.
   * 
   * Space is allocated as for a diagonal matrix and elements are set to 1.
   * This storage could be saved.
   */ 

class IdentityMatrix : public DiagonalMatrix
{

public :
   
    // = CONSTRUCTORS

    IdentityMatrix (int n) ;

    // = INQUIRY

    inline MatrixType GetType () ;
} ;


/*! \class TensorMatrix
   * \brief Manages tensor matrices.
   * 
   * Tensor matrices are made of directional matrices. The full matrix
   * is never provided nor built explicitly. Consequently, this matrix
   * class encapsulates its directional matrices instead of deriving from
   * the Matrix base class.
   * The directional matrices are <{matrix1}>, <{matrix2}> and <{matrix3}>.
   * They are stored separately, rather than in a Vector<Matrix*> object, 
   * for efficiency reasons: this avoids costly allocations/deallocations
   * of the vector itself. Also for efficiency reasons, when a tensor
   * matrix is assigned a directional matrix m, it simply points to m
   * instead of making a copy of m.
   * 
   * The dimension <{dim}>, or equivalently the number of directional
   * matrices, is limited to 3.
   * 
   * Operations of multiplication by a vector are provided: in 3D, these
   * are expressed by:
   *   y = (A3 o A2 o A1) x   or  y = (A3(T) o A2(T) O A1(T)) x ,
   * where A_i stands for <{matrix_i}>, "o" for the tensor product and "T"
   * for a matrix tranposition.
   * 
   * Since tensor matrices are frequently allocated/deallocated, pseudo-
   * dynamic allocations/deallocations of tensor matrices are handled by a
   * <{pool}>.
   * 
   * Furthermore, profiling tests show that dynamic allocations/
   * deallocations of temporary full matrices in methods <{MatVec}> and
   * <{MatTransVec}> are not cost effective. Therefore two matrices
   * <{auxiliaryMatrix1}> and <{auxiliaryMatrix1}> have been defined as
   * static attribute of the class.
   */ 
   
class TensorMatrix
{ 

protected :

    int     dimension ;
    Matrix* matrix1 ;
    Matrix* matrix2 ;
    Matrix* matrix3 ;

    static Pool*       pool ;
    static FullMatrix* auxiliaryMatrix1 ;
    static FullMatrix* auxiliaryMatrix2 ;

public :


    // = PSEUDO-DYNAMIC ALLOCATION

    inline void* operator new (size_t) ;
    inline void  operator delete (void* A) ;

    // = CONSTRUCTORS

    inline TensorMatrix (int dim) ;

    // ACCESS TO THE MATRICES

    void            SetMatrix (int i, Matrix* mat) ;
    void            DeleteMatrices () ;

    // = INQUIRY

    Matrix*         GetMatrix (int i) ;
    int             GetNrow () ;
    int             GetNcol () ;
    boolean         IsConsistent () ;
    void            Print () ;

    // = OPERATIONS

    DiagonalMatrix* Inverse (boolean modifyZeroEigenvalues=false) ;
    void            MatVec (RealVector* x, RealVector* y) ;
    void            MatTransVec (RealVector* x, RealVector* y) ;

    // = DEBUGGING

    FullMatrix*     Expanded () ;
} ;



//_______________________________ Matrix (inline) _____________________________


Matrix :: Matrix (int n, int m) 
{
    // Constructor. Initializes the receiver to a matrix of size (n,m).

# ifdef REQUIRE
    Require("valid number of rows 'n'", n >= 0) ;
    Require("valid number of columns 'm'", m >= 0) ;
    InMethod("Matrix::Matrix(n,m)") ;
# endif

    nrow = n ;
    ncol = m ;
}


Matrix :: ~Matrix ()
{
    // Destructor.
}


int Matrix :: GetNcol ()
{
    // Returns the number of columns of the receiver.

    return ncol ; 
}


int Matrix :: GetNrow ()
{
    // Returns the number of rows of the receiver.

    return nrow ; 
}


real* Matrix :: GetValues ()
{
    // Returns a pointer to the first coefficient of the receiver.

    return values ; 
}


//_____________________________ FullMatrix (inline) ___________________________


FullMatrix :: FullMatrix () 
    : Matrix (0,0)
{
    // Constructor. Initializes the receiver to a full matrix of size (0,0). 

    values        = NULL ;
    allocatedSize = 0 ;
    isAlias       = false ;
}   


FullMatrix :: FullMatrix (int n, int m)
    : Matrix (n,m)
{
    // Constructor. Initializes the receiver to a full matrix of size (0,0).

    allocatedSize = n*m ;
    if (allocatedSize)
	values = new real[allocatedSize] ;
    else
	values = NULL ;

    isAlias = false ;
}


FullMatrix :: FullMatrix (RealVector* x, int n, int m) 
    : Matrix (n,m)
{
    // Constructor. Initializes the receiver as an alias of size (n,m) to 'x',
    // i.e., uses the values of 'x' as values of the receiver.

# ifdef REQUIRE
    Require("'x' is large enough", x->GetSize() >= n*m) ;
    InMethod("FullMatrix(x,n,m)") ;
# endif

    values        = NULL ;
    isAlias       = true ;
    DefineAlias(x,n,m) ;
    allocatedSize = x->GetSize() ;
}


FullMatrix :: ~FullMatrix ()
{
    // Destructor.
  
      if (! isAlias)
      	delete values ;
}


void FullMatrix :: DefineAlias (RealVector* x, int n, int m)
{
    // Aliases the receiver matrix to 'x', with 'n' rows and 'm' columns.
    // The values of 'x' are used as values of the receiver.

# ifdef REQUIRE
    Require("size(A) < size(x)", n*m <= x->GetSize()) ;
    InMethod("FullMatrix::DefineAlias(x,n,m)") ;
# endif

    values        = x->GetValues() ;
    nrow          = n ;
    ncol          = m ;
    allocatedSize = x->GetSize() ;
    isAlias       = true ;
}


int FullMatrix :: GetAllocatedSize ()
{
    // Returns the number of coefficients physically allocated.
    // This number is >= nrow*ncol.

    return allocatedSize ;
}


real* FullMatrix :: GetColumn (int j)
{
    // Returns a pointer to the beginning of the 'j'-th column of the receiver.

# ifdef REQUIRE
    Require("valid argument 'j'", j>=1 && j<=ncol) ;
    InMethod("FullMatrix::GetColumn(j)") ;
# endif

    return &values[(j-1)*nrow] ; 
}


MatrixType FullMatrix :: GetType ()
{
    // Returns the type of the receiver.
 
    return FULL_MT ; 
}

void FullMatrix :: SetDimensions (int n, int m)
{
    // Resets the dimensions of the receiver to n*m.

# ifdef REQUIRE
    Require("Valid new dimension", n*m >= 1 && n*m <= allocatedSize) ;
    InMethod("FullMatrix::SetDimensions(n,m)") ;
# endif

    nrow = n ;
    ncol = m ;
}


//_____________________________ DiagonalMatrix (inline) _______________________


MatrixType DiagonalMatrix :: GetType () 
{
    // Returns the type of the receiver.

    return DIAGONAL_MT ; 
}


//____________________________ IdentityMatrix (inline) ________________________


MatrixType IdentityMatrix :: GetType ()
{
    // Returns the type of the receiver.
 
    return IDENTITY_MT ; 
}


//___________________________ TensorMatrix (inline) ___________________________


void* TensorMatrix :: operator new (size_t)
{
    // Allocation operator.

    return pool->Alloc() ;
}


void TensorMatrix :: operator delete (void* A)
{
    // Deallocation operator.

    pool->Free(A) ;
}


TensorMatrix :: TensorMatrix (int dim)
{
    // Constructor. Initializes the receiver to a tensor matrix of order 'dim'.
    // For efficiency purposes, the 'dim' directional matrices are not allocated.

# ifdef REQUIRE
    Require("Valid order 'dim'", dim >= 1 && dim <= 3) ;
    InMethod("TensorMatrix::TensorMatrix(dim)") ;
# endif

    dimension = dim ;
    matrix1   = NULL ;
    matrix2   = NULL ;
    matrix3   = NULL ;
}


#endif

