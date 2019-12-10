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

#ifndef PARAM
#define PARAM

/*!
 * \file param.hxx
 * \brief This file contains some parameters
 * \author {Lots of authors}
 * \version 1.0
 */


using namespace std ;
#include <iostream>


// ------------
// real NUMBERS
// ------------

typedef double real ;

const real ZERO               =  0. ;
const real ONE                =  1. ;
const real TWO                =  2. ;
const real THREE              =  3. ;
const real FOUR               =  4. ;
const real FIVE               =  5. ;
const real SIX                =  6. ;
const real SEVEN              =  7. ;
const real EIGHT              =  8. ;
const real ONE_MINUS          = -1. ;
const real P1                 =  ONE ; 
const real M1                 =  ONE_MINUS ;
const real HALF               =  0.5 ;
const real THIRD              =  ONE / THREE ;
const real ONE_SEVENTH        =  ONE / SEVEN ;
const real EIGHT_SEVENTH      =  ONE / EIGHT ;
const real TWO_THIRD          =  TWO / THREE ;
const real FIVE_HALF          =  FIVE / HALF ;
const real SIXTH              =  ONE / SIX ;
const real DEFAULT_TOLERANCE  =  1.e-12 ;
const real UNINITIALIZED_REAL =  123456789. ;
const real PI                 =  3.14159265358979323846264338327950 ;
//  avoids including math.h!

// --------
// INTEGERS
// --------

const int IZERO                 =  0 ;
const int IONE                  =  1 ;
const int UNINITIALIZED_INTEGER = -99999999 ;
const int ALL                   =  0 ;
const int MAXWORD               = 256 ;

// --------
// BOOLEANS
// --------

typedef int boolean ;

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif


// ----------------
// THE NULL POINTER
// ----------------

#ifndef NULL
#define NULL 0
#endif


// --------------------
// ADDITIONAL FUNCTIONS
// --------------------

#ifndef maxfct
#define maxfct(a,b) ( a > b ? a : b )
#endif

#ifndef minfct
#define minfct(a,b) ( a < b ? a : b )
#endif


// -----------------------------------
// FORTRAN INTERFACE (BLAS AND LAPACK)
// -----------------------------------

#if defined SILICON_GRAPHICS
#  define DAXPY  daxpy_
#  define DCOPY  dcopy_
#  define DDOT   ddot_
#  define DNRM2  dnrm2_
#  define DSCAL  dscal_
#  define DGEMV  dgemv_
#  define DGEMM  dgemm_
#  define DGEEV  dgeev_
#  define DGESV  dgesv_
#  define DGETRF dgetrf_
#  define DGETRI dgetri_
#  define DGETRS dgetrs_
#  define DGEESX dgeesx_
#elif defined HEWLETT_PACKARD
#  define DAXPY  daxpy
#  define DCOPY  dcopy
#  define DDOT   ddot
#  define DNRM2  dnrm2
#  define DSCAL  dscal
#  define DGEMV  dgemv
#  define DGEMM  dgemm
#  define DGEEV  dgeev
#  define DGESV  dgesv
#  define DGETRF dgetrf
#  define DGETRI dgetri
#  define DGETRS dgetrs
#  define DGEESX dgeesx
#elif defined DIGITAL
#  define DAXPY  daxpy_
#  define DCOPY  dcopy_
#  define DDOT   ddot_
#  define DNRM2  dnrm2_
#  define DSCAL  dscal_
#  define DGEMV  dgemv_
#  define DGEMM  dgemm_
#  define DGEEV  dgeev_
#  define DGESV  dgesv_
#  define DGETRF dgetrf_
#  define DGETRI dgetri_
#  define DGETRS dgetrs_
#  define DGEESX dgeesx_
#elif defined SUN
#  define DAXPY  daxpy_
#  define DCOPY  dcopy_
#  define DDOT   ddot_
#  define DNRM2  dnrm2_
#  define DSCAL  dscal_
#  define DGEMV  dgemv_
#  define DGEMM  dgemm_
#  define DGEEV  dgeev_
#  define DGESV  dgesv_
#  define DGETRF dgetrf_
#  define DGETRI dgetri_
#  define DGETRS dgetrs_
#  define DGEESX dgeesx_
#elif defined CRAY_T3D
#  define DAXPY  SAXPY
#  define DCOPY  SCOPY
#  define DDOT   SDOT
#  define DNRM2  SNRM2
#  define DSCAL  SSCAL
#  define DGEMV  SGEMV
#  define DGEMM  SGEMM
#  define DGEEV  SGEEV
#  define DGESV  SGESV
#  define DGETRF SGETRF
#  define DGETRI SGETRI
#  define DGETRS SGETRS
#  define DGEESX SGEESX
#elif defined GNU
#  define DAXPY  daxpy_
#  define DCOPY  dcopy_
#  define DDOT   ddot_
#  define DNRM2  dnrm2_
#  define DSCAL  dscal_
#  define DGEMV  dgemv_
#  define DGEMM  dgemm_
#  define DGEEV  dgeev_
#  define DGESV  dgesv_
#  define DGETRF dgetrf_
#  define DGETRI dgetri_
#  define DGETRS dgetrs_
#  define DGEESX dgeesx_
#elif defined BLUE_GENE
#  define DAXPY  daxpy
#  define DCOPY  dcopy
#  define DDOT   ddot
#  define DNRM2  dnrm2
#  define DSCAL  dscal
#  define DGEMV  dgemv
#  define DGEMM  dgemm
#  define DGEEV  dgeev
#  define DGESV  dgesv
#  define DGETRF dgetrf
#  define DGETRI dgetri
#  define DGETRS dgetrs
#  define DGEESX dgeesx
#endif


extern "C" {

    // BLAS 1
    void DAXPY  (int* n, real* alpha, real* x, int* incx, real* y, int* incy) ;
    void DCOPY  (int* n, real* x, int* incx, real* y, int* incy) ;
    real DDOT   (int* n, real* x, int* incx, real* y, int* incy) ;
    real DNRM2  (int* n, real* x, int* incx) ;
    void DSCAL  (int* n, real* alpha, real* x, int* incx) ;

    // BLAS 2
    void DGEMV  (char* trans, int* m, int* n, real* alpha, real* a, int* lda,real* x, int* incx, real* beta, real* y, int* incy) ;
    // BLAS 3
    void DGEMM  (char* transa, char* transb, int* m, int* n, int* k, real* alpha,real* a, int* lda, real* b, int* ldb, real* beta,real* c, int* ldc) ;
    // LAPACK
    void DGEEV  (char* jobvl, char* jobvr, int* n, real* a, int* lda,real* wr, real* wi, real* vl, int* ldvl, real* vr, int* ldvr,
		 real* work, int* lwork, int* info) ;
    void DGESV  (int* n, int* nrhs, real* a, int* lda, int* ipiv,real* b, int* ldb, int* info) ;
    void DGETRF (int* m, int* n, real* a, int* lda, int* ipiv, int* info) ;
    void DGETRI (int* n, real* a, int* lda, int* ipiv, real* work, int* lwork,int* info) ;
    void DGETRS (char* trans, int* n, int* nrhs, real* a, int* lda, int* ipiv,real* b, int* ldb, int* info) ;

    void DGEESX(char* jobsc, char* SORT, boolean* SELECT, char* SENSE,int* N, real *A,int* LDA, int* SDIM,real* WR, real* WI, real* VS,int* LDVS,real* RCONDE,real* RCONDV, real* WORK,int* LWORK,int* IWORK,int* LIWORK,boolean* BWORK,int* INFO ) ;

}

#endif
