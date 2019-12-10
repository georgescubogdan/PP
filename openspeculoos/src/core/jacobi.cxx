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

// jacobi.cxx

#include "core/jacobi.hxx"
#include "core/matrix.hxx"
#include "core/vector.hxx"
#include "core/vector_templates.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include <math.h>


void DGJD (FullMatrix* D, RealVector* Z, real ALPHA, real BETA)
{
    // Computes the derivative matrix D associated with the Nth order Lagrangian
    // interpolants through the NZ Gauss Jacobi points contained in Z.

# ifdef REQUIRE
    Require("the matrix is squared", D->GetNcol() == D->GetNrow()) ;
    Require("matrix order = vector size", D->GetNcol() == Z->GetSize()) ;
    Require("ALPHA > -1", ALPHA > -1) ;
    Require("BETA > -1", BETA > -1) ;
    InMethod("DGJD(D,Z,ALPHA,BETA)") ;
# endif

# ifdef CHECK
    if (Z->GetSize() == 1)
	Warning("function DGJD",
		"case Gauss-Legendre of degree 0: returns the nul matrix") ;
# endif

    real PPI,PDI,PM1,PDM1,PM2,PDM2,PPJ,PDJ ;
    int  I,J,NZ ;

    NZ = Z->GetSize() ;
    for (J=1 ; J<=NZ ; J++) {
	for (I=1 ; I<=NZ ; I++) {
	    JACOBF(PPI,PDI,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z->At(I)) ;
	    JACOBF(PPJ,PDJ,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z->At(J)) ;
	    if (I!=J) 
		D->At(I,J) = PDI / (PDJ*(Z->At(I)-Z->At(J))) ;
	    if (I==J) 
		D->At(I,J) = ((ALPHA+BETA+TWO)*Z->At(I)+ALPHA-BETA) /
		    (TWO*(ONE-Z->At(I)*Z->At(I))) ;
	}
    }
}


void DGLJD (FullMatrix* D, RealVector* Z, real ALPHA, real BETA)
{
    // Computes the derivative matrix D associated with the Nth order 
    // Lagrangian interpolants through the NZ Gauss-Lobatto Jacobi points 
    // contained in Z.
    // Note: D must be a square matrix and its order must be equal to the size
    //       of Z

# ifdef REQUIRE
    Require("D is squared", D->GetNcol() == D->GetNrow()) ;
    Require("order of D = size of Z", D->GetNcol() == Z->GetSize()) ;
    Require("minimum number of Gauss-Lobatto points is 2", Z->GetSize() >= 2) ;
    Require("ALPHA > -1", ALPHA > -1) ;
    Require("BETA > -1", BETA > -1) ;
    InMethod("DGLJD(D,Z,ALPHA,BETA)") ;
#endif

    int I,J ;

    int NZ      = Z->GetSize() ;
    int N       = NZ - 1 ;
    real DN     = (real) N ;
    real EIGVAL = -DN*(DN+ALPHA+BETA+ONE) ;

    real CI,CJ,PPI,PDI,PM1,PDM1,PM2,PDM2,PPJ,PDJ ;

    for (J=1 ; J<=NZ ; J++) {
	for (I=1 ; I<=NZ ; I++) {
	    JACOBF(PPI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z->At(I)) ;
	    JACOBF(PPJ,PDJ,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z->At(J)) ;
	    CI = EIGVAL*PPI-(BETA*(ONE-Z->At(I))-ALPHA*(ONE+Z->At(I)))*PDI ;
	    CJ = EIGVAL*PPJ-(BETA*(ONE-Z->At(J))-ALPHA*(ONE+Z->At(J)))*PDJ ;
	    if (I != J)  
		D->At(I,J) = CI/(CJ*(Z->At(I)-Z->At(J))) ;
	    else {
		if ((I != 1) && (I != NZ))
		    D->At(I,J) = (ALPHA*(ONE+Z->At(I))-BETA*(ONE-Z->At(I))) /
			(TWO*(ONE-Z->At(I)*Z->At(I))) ;
		else if (I == 1)
		    D->At(I,J) =  (EIGVAL+ALPHA) / (TWO*(BETA+TWO)) ;
		else if (I == NZ)
		    D->At(I,J) = -(EIGVAL+BETA) / (TWO*(ALPHA+TWO)) ;
	    }
	}
    }
}


void IGJMD (FullMatrix* I12, RealVector* Z1, RealVector* Z2, 
	    real ALPHA, real BETA)
{
    // Computes the one-dimensional interpolation operator (matrix) I12
    // for interpolating a variable from a Gauss Jacobi mesh (1)  
    // to a another mesh M (2).
    // Z1 : NZ1 Gauss Jacobi points.
    // Z2 : NZ2 points on mesh M.

# ifdef REQUIRE
    Require("nrow(I12) = size(Z2)", I12->GetNrow() == Z2->GetSize()) ;
    Require("ncol(I12) = size(Z1)", I12->GetNcol() == Z1->GetSize()) ;
    Require("ALPHA > -1", ALPHA > -1) ;
    Require("BETA > -1", BETA > -1) ;
    InMethod("IGJMD(I12,Z1,Z2,ALPHA,BETA)") ;
# endif

    int NZ1 = Z1->GetSize() ;
    int NZ2 = Z2->GetSize() ;

    real ZI, ZJ, DZ ;
    real EPS = 1.e-5;

    int I,J ;
    real PPI,PDI,PM1,PDM1,PM2,PDM2,PPJ,PDJ ;

    if (NZ1 == 1)
	// YDP modif: error in the next line in all previous versions of VVK and 
	// YDP: there used to be 'I12(1,1) = 1. ;' instead of the following loop on
	// J.
	for (J=1 ; J<=NZ2 ; J++)
	    I12->At(J,1) = 1. ;

    else {
	for (J=1 ; J<=NZ1 ; J++) {
	    ZJ = Z1->At(J) ;
	    for (I=1 ; I<=NZ2 ; I++) {
		ZI = Z2->At(I) ;
		DZ = ZI - ZJ ;
		if (fabs(DZ) < EPS) 
		    I12->At(I,J) = 1. ;
		else {
		    JACOBF(PPI,PDI,PM1,PDM1,PM2,PDM2,NZ1,ALPHA,BETA,ZI) ;
		    JACOBF(PPJ,PDJ,PM1,PDM1,PM2,PDM2,NZ1,ALPHA,BETA,ZJ) ;
		    I12->At(I,J) = PPI / (PDJ*(ZI-ZJ)) ;
		}
	    }
	}
    }
}


void IGLJMD (FullMatrix* I12, RealVector* Z1, RealVector* Z2, 
	     real ALPHA, real BETA, int interior)
{
    // Computes the one-dimensional interpolation operator (matrix) I12
    // for interpolating a variable from a Gauss-Lobatto Jacobi mesh (1)  
    // to a another mesh M (2).
    // Z1 : NZ1 Gauss Jacobi points.
    // Z2 : NZ2 points on mesh M.
    // If interior != 0 then only the interior collocation points are considered
    // i.e. all except {-1,1}.

# ifdef REQUIRE
    Require("nrow(I12) = size(Z2)", I12->GetNrow() == Z2->GetSize()) ;
    Require("ncol(I12) = size(Z1)", I12->GetNcol() == Z1->GetSize()) ;
    Require("ALPHA > -1", ALPHA > -1) ;
    Require("BETA > -1", BETA > -1) ;
    InMethod("IGLJMD(I12,Z1,Z2,ALPHA,BETA,interior)") ;
# endif

    int NZ1 = Z1->GetSize() ;
    int NZ2 = Z2->GetSize() ;

    real ZI, ZJ, DZ ;
    real EPS = 1.e-5;

    int N ; 
    real DN, EIGVAL, fac1, fac2 ;

    if (interior) 
	N = NZ1 + 1 ;
    else
	N = NZ1 - 1 ;

    DN = (real) N ;
    EIGVAL = -DN*(DN+ALPHA+BETA+1.) ;

    int I,J ;
    real PPI,PDI,PM1,PDM1,PM2,PDM2,PPJ,PDJ ;

    for (J=1 ; J<=NZ1 ; J++) {
	ZJ = Z1->At(J) ;
	for (I=1 ; I<=NZ2 ; I++) {
	    ZI = Z2->At(I) ;
	    DZ = ZI - ZJ ;
	    if (fabs(DZ) < EPS) 
		I12->At(I,J) = 1. ;
	    else {
		if (interior) fac1 = (1. - ZJ*ZJ) ;
		else          fac1 = (1. - ZI*ZI) ;
		JACOBF(PPI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZI) ;
		JACOBF(PPJ,PDJ,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZJ) ;
		fac2 = EIGVAL*PPJ+ALPHA*(1.+ZJ)*PDJ-BETA*(1.-ZJ)*PDJ ;
		I12->At(I,J) = fac1*PDI / (fac2*(ZI-ZJ)) ;
	    }
	}
    }
}


void PHIML (FullMatrix* mat, RealVector* z)
{
    // MAH 30.11.07
    
    real zi ;
    int  i, j, n, size ;
      
    size = z->GetSize() ;      
    
    for (i=1 ; i<=size ; i++) {
	zi = z->At(i) ;

	mat->At(i,1) = (1 - zi) / 2. ;                                        // Phi_0(ksi_i)
	mat->At(i,2) = (1 + zi) / 2. ;                                        // Phi_1(ksi_i)

	for (j=3 ; j<=size ; j++) mat->At(i,j) = PNLEG(zi,j-1)-PNLEG(zi,j-3) ;// Phi_j(ksi_i) = L_j(ksi_i) - L_j-2(ksi)
    }
}


real PNLEG (real z, int n)
{
    // MAH 30.11.07

    int  k ;
    real p1, p2, p3, pnleg, fk ;
    
    if (fabs(z) < 1.0e-25) z = 0.0 ;
    
    p1 = 1. ;
    if (n==0) { 
	pnleg = p1 ;
	return pnleg ;
         }

    p2 = z ;
    if (n==1) { 
	pnleg = p2 ;
	return pnleg ;
         }

     
    for (k=1 ; k<=n-1 ; k++) {
	fk = real(k) ;
	p3 = ((2.*fk+1.)*z*p2 - fk*p1)/(fk+1.) ;
	p1 = p2 ;
	p2 = p3 ;
    }
    pnleg = p3 ;
    return pnleg ;
}


void ZIGLJD (RealVector* Z, real ALPHA, real BETA)
{
    // Generates NP=size(Z) GAUSS-LOBATTO JACOBI points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP+1 i.e. coordinates are given inside 
    // the interval ]-1,1[.

    int NP = Z->GetSize() ;

    real ALPG,BETG,APB ;

    ALPG  = ALPHA+1. ;
    BETG  = BETA+1. ;
    APB   = ALPG+BETG ;
  
    if (NP == 1)
	Z->At(1) = (BETA-ALPHA)/(APB+2.) ;
    else
	JACG(Z,ALPG,BETG) ;
}


void ZWGLJD (RealVector* Z, RealVector* W, real ALPHA, real BETA)
{
    // Generate NP=size(Z) GAUSS-LOBATTO JACOBI points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP-1.

    int I,NP = Z->GetSize() ;

    int N   = NP-1 ;
    int NM1 = N-1 ;

    real ALPG,BETG ;

    RealVector workZ,workW ;

    real P,PD,PM1,PDM1,PM2,PDM2 ;

    if (NM1 > 0) {
	ALPG  = ALPHA+1. ;
	BETG  = BETA+1. ;
	workZ.DefineAlias(Z,2,NM1) ;
	workW.DefineAlias(W,2,NM1) ;
	ZWGJD(&workZ,&workW,ALPG,BETG) ;
    }

    Z->At(1)  = -1. ;
    Z->At(NP) = 1. ;

    for (I=2 ; I<=NP-1 ; I++) 
	W->At(I) = W->At(I) / (1.-Z->At(I)*Z->At(I)) ;

    if (NM1 > 0) {
	JACOBF(P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z->At(1)) ;
	W->At(1)   = ENDW1(N,ALPHA,BETA)/(2.*PD) ;
	JACOBF(P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z->At(NP)) ;
	W->At(NP)  = ENDW2(N,ALPHA,BETA)/(2.*PD) ;
    } else {
	W->At(1)  = 1. ;
	W->At(NP) = 1. ;
    }
}


void ZGLJD (RealVector* Z, real ALPHA, real BETA)
{
    // Generate NP=size(Z) GAUSS-LOBATTO JACOBI points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP-1.

    int I,NP = Z->GetSize() ;

    int N   = NP-1 ;
    int NM1 = N-1 ;

    real ALPG,BETG ;

    RealVector workZ ;

    real P,PD,PM1,PDM1,PM2,PDM2 ;

    if (NM1 > 0) {
	ALPG  = ALPHA+1. ;
	BETG  = BETA+1. ;
	workZ.DefineAlias(Z,2,NM1) ;
	ZGJD(&workZ,ALPG,BETG) ;
    }

    Z->At(1)  = -1. ;
    Z->At(NP) = 1. ;
}


void ZWIGLJD (RealVector* Z, RealVector* ZW, RealVector* W, real ALPHA, real BETA)
{
    // Generates NP=size(Z) GAUSS-LOBATTO JACOBI inner points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP+1 i.e. coordinates are given inside 
    // the interval ]-1,1[.

    int NP = Z->GetSize() ;

    real ALPG,BETG,APB ;

    ALPG  = ALPHA+1. ;
    BETG  = BETA+1. ;
    APB   = ALPG+BETG ;
  
    if (NP == 1)
	Z->At(1) = (BETA-ALPHA)/(APB+2.) ;
    else
	JACG(Z,ALPG,BETG) ;  
    
    // Generates NPW=size(Z)+2 GAUSS-LOBATTO JACOBI weights (W) and points
    // (ZW) associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).

    int I,NPW = ZW->GetSize();  

    int N   = NPW-1 ;
    int NM1 = N-1 ;

    RealVector workZ,workW ;

    real P,PD,PM1,PDM1,PM2,PDM2 ;

    if (NM1 > 0) {
	ALPG  = ALPHA+1. ;
	BETG  = BETA+1. ;
	workZ.DefineAlias(ZW,2,NM1) ;
	workW.DefineAlias(W,2,NM1) ;
	ZWGJD(&workZ,&workW,ALPG,BETG) ;
    }

    ZW->At(1)  = -1. ;
    ZW->At(NPW) = 1. ;

    for (I=2 ; I<=NPW-1 ; I++) 
	W->At(I) = W->At(I) / (1.-ZW->At(I)*ZW->At(I)) ;

    if (NM1 > 0) {
	JACOBF(P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZW->At(1)) ;
	W->At(1)   = ENDW1(N,ALPHA,BETA)/(2.*PD) ;
	JACOBF(P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZW->At(NPW)) ;
	W->At(NPW)  = ENDW2(N,ALPHA,BETA)/(2.*PD) ;
    } else {
	W->At(1)   = 1. ;
	W->At(NPW) = 1. ;
    }
  
}


void ZWGJD (RealVector* Z, RealVector* W, real ALPHA, real BETA)
{
    // Generate NP=size(Z) GAUSS JACOBI points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP-1.

    int I,NP = Z->GetSize() ;

    int N   = NP-1 ;

    real APB = ALPHA+BETA ;

    RealVector workZ,workW ;

    real P,PD,PM1,PDM1,PM2,PDM2 ;

    if (NP == 1) {

	Z->At(1) = (BETA-ALPHA)/(APB+2.) ;
	W->At(1) = GAMMAF(ALPHA+1.)*GAMMAF(BETA+1.)/GAMMAF(APB+2.)*pow(2.,APB+1.) ;

    } else {

	JACG(Z,ALPHA,BETA) ;
    
	int NP1   = N+1 ;
	int NP2   = N+2 ;
	real DNP1  = (real) NP1 ;
	real DNP2  = (real) NP2 ;
	real FAC1  = DNP1+ALPHA+BETA+1. ;
	real FAC2  = FAC1+DNP1 ;
	real FAC3  = FAC2+1. ;
	real FNORM = PNORMJ(NP1,ALPHA,BETA) ;
	real RCOEF = (FNORM*FAC2*FAC3)/(2.*FAC1*DNP2) ;

	for (I=1 ; I<=NP ; I++) {
	    JACOBF(P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z->At(I)) ;
	    W->At(I)  = -RCOEF/(P*PDM1) ;
	}
    }
}

void ZGJD (RealVector* Z, real ALPHA, real BETA)
{
    // Generate NP=size(Z) GAUSS JACOBI points (Z)
    // associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
    // The polynomial degree N=NP-1.

    int I,NP = Z->GetSize() ;

    int N   = NP-1 ;

    real APB = ALPHA+BETA ;

    RealVector workZ ;

    real P,PD,PM1,PDM1,PM2,PDM2 ;

    if (NP == 1) {

	Z->At(1) = (BETA-ALPHA)/(APB+2.) ;

    } else {

	JACG(Z,ALPHA,BETA) ;
    
// 	int NP1   = N+1 ;
// 	int NP2   = N+2 ;
// 	real DNP1  = (real) NP1 ;
// 	real DNP2  = (real) NP2 ;
// 	real FAC1  = DNP1+ALPHA+BETA+1. ;
// 	real FAC2  = FAC1+DNP1 ;
// 	real FAC3  = FAC2+1. ;
// 	real FNORM = PNORMJ(NP1,ALPHA,BETA) ;
// 	real RCOEF = (FNORM*FAC2*FAC3)/(2.*FAC1*DNP2) ;
    }
}


// ------------------
// Internal functions
// ------------------

real ENDW1 (int N, real ALPHA, real BETA)
{
    real APB   = ALPHA+BETA ;

    real A1,A2,A3,F1,F2,F3,FINT1,FINT2,retval ;
    real DI,ABN,ABNN ;

    if (N == 0) 
	retval = ZERO ;
    else {
	F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE) ;
	F1   = F1*(APB+TWO)*pow(2.,APB+TWO)/TWO ;
	if (N == 1) 
	    retval = F1 ;
	else {
	    FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE) ;
	    FINT1 = FINT1*pow(2.,APB+TWO) ;
	    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR) ;
	    FINT2 = FINT2*pow(2.,APB+THREE) ;
	    F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2)* (APB+THREE)/FOUR ;
	    if (N == 2) 
		retval = F2 ;
	    else {
		for (int I=3 ; I<=N ; I++) {
		    DI   = (real) (I-1) ;
		    ABN  = ALPHA+BETA+DI ;
		    ABNN = ABN+DI ;
		    A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE)) ;
		    A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO)) ;
		    A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE)) ;
		    F3   =  -(A2*F2+A1*F1)/A3 ;
		    F1   = F2 ;
		    F2   = F3 ;
		}
		retval  = F3 ;
	    }
	}
    }

    return retval ;
}


real ENDW2 (int N, real ALPHA, real BETA)
{
    real APB   = ALPHA+BETA ;

    real A1,A2,A3,F1,F2,F3,FINT1,FINT2,retval ;
    real DI,ABN,ABNN ;

    APB   = ALPHA+BETA ;
    if (N == 0)
	retval = ZERO ;
    else {
	F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE) ;
	F1   = F1*(APB+TWO)*pow(2.,APB+TWO)/TWO ;
	if (N == 1)
	    retval = F1 ;
	else {
	    FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE) ;
	    FINT1 = FINT1*pow(2.,APB+TWO) ;
	    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR) ;
	    FINT2 = FINT2*pow(2.,APB+THREE) ;
	    F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2) * (APB+THREE)/FOUR ;
	    if (N == 2)
		retval = F2 ;
	    else {
		for (int I=3; I<=N ; I++) {
		    DI   = (real) (I-1) ;
		    ABN  = ALPHA+BETA+DI ;
		    ABNN = ABN+DI ;
		    A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE)) ;
		    A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO)) ;
		    A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE)) ;
		    F3   =  -(A2*F2+A1*F1)/A3 ;
		    F1   = F2 ;
		    F2   = F3 ;
		}
		retval  = F3 ;
	    }
	}
    }
    return retval ;
}


real GAMMAF(real X)
{
    real retval = 1. ;

    if (X ==-0.5) retval =  -2.*sqrt(PI) ;
    if (X == 0.5) retval =   sqrt(PI) ;
    if (X == 1. ) retval =   1. ;
    if (X == 2. ) retval =   1. ;
    if (X == 1.5) retval =   sqrt(PI) / 2. ;
    if (X == 2.5) retval =   1.5 * sqrt(PI) / 2. ;
    if (X == 3.5) retval =   2.5 * 1.5 * sqrt(PI) / 2. ;
    if (X == 3. ) retval =   2. ;
    if (X == 4. ) retval =   6. ;
    if (X == 5. ) retval =  24. ;
    if (X == 6. ) retval = 120. ;

    return retval ;
}


void JACG (RealVector* XJAC, real ALPHA, real BETA)
{
    // Compute NP Gauss points XJAC, which are the zeros of the
    // Jacobi polynomial J(NP) with parameters ALPHA and BETA.
    // ALPHA and BETA determines the specific type of Gauss points.
    // Examples:
    //  ALPHA = BETA =  0.0  ->  Legendre points
    //  ALPHA = BETA = -0.5  ->  Chebyshev points

    int I,J,K,JM,JMIN ;
    real P,PD,PM1,PDM1,PM2,PDM2,X,X1,X2,DELX,XLAST,RECSUM,XMIN,SWAP ;

    int KSTOP = 10 ;
    real EPS = 1.0e-12 ;

    int  NP  = XJAC->GetSize() ;
    int  N   = NP-1 ;
    real DTH = PI/((real) (2*N+2)) ;

    XLAST = ZERO ;

    for (J=1 ; J<=NP ; J++) {
	if (J == 1) 
	    X  = cos(((real)(2*(J-1)+1)) * DTH) ;
	else {
	    X1 = cos(((real)(2*(J-1)+1)) * DTH) ;
	    X2 = XLAST ;
	    X  = (X1+X2)/2. ;
	}
	for (K=1 ; K<=KSTOP ; K++) {
	    JACOBF(P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X) ;
	    RECSUM = 0. ;
	    JM = J-1 ;
	    for (I=1 ; I<=JM ; I++)
		RECSUM = RECSUM+1./(X-XJAC->At(NP-I+1)) ;
	    DELX = -P/(PD-RECSUM*P) ;
	    X    = X+DELX ;
	    if (fabs(DELX) < EPS) break ;
	}
	XJAC->At(NP-J+1) = X ;
	XLAST        = X ;
    }
    for (I=1 ; I<=NP ; I++) {
	XMIN = 2. ;
	for (J=I ; J<=NP ; J++)
	    if (XJAC->At(J) < XMIN) {
		XMIN = XJAC->At(J) ;
		JMIN = J ;
	    }
	if (JMIN != I) {
	    SWAP = XJAC->At(I) ;
	    XJAC->At(I) = XJAC->At(JMIN) ;
	    XJAC->At(JMIN) = SWAP ;
	}
    }
}


void JACOBF(real &POLY, real &PDER, real &POLYM1, real &PDERM1, real &POLYM2,
	    real &PDERM2, int N, real ALP, real BET, real X)
{
    // Computes the Jacobi polynomial (POLY) and its derivative (PDER)
    // of degree N at X.

    real DK,A1,A2,A3,A4,B3,POLYN,PDERN,PSAVE,PDSAVE ;

    real APB  = ALP+BET ;

    POLY = 1. ;
    PDER = 0. ;

    if (N != 0) {
	real POLYL = POLY ;
	real PDERL = PDER ;
    
	POLY  = (ALP-BET+(APB+2.)*X)/2. ;
	PDER  = (APB+2.)/2. ;
	if (N != 1) {
	    for (int K=2 ; K<=N ; K++) {
		DK = (real) K ;
		A1 = 2.*DK*(DK+APB)*(2.*DK+APB-2.) ;
		A2 = (2.*DK+APB-1.)*(ALP*ALP-BET*BET) ;
		B3 = (2.*DK+APB-2.) ;
		A3 = B3*(B3+1.)*(B3+2.) ;
		A4 = 2.*(DK+ALP-1.)*(DK+BET-1.)*(2.*DK+APB) ;
		POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1 ;
		PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1 ;
		PSAVE  = POLYL ;
		PDSAVE = PDERL ;
		POLYL  = POLY ;
		POLY   = POLYN ;
		PDERL  = PDER ;
		PDER   = PDERN ;
	    }
	    POLYM1 = POLYL ;
	    PDERM1 = PDERL ;
	    POLYM2 = PSAVE ;
	    PDERM2 = PDSAVE ;
	}
    }
}


real PNORMJ(int N, real ALPHA, real BETA)
{
    real DN    = (real) N ;
    real CONST = ALPHA+BETA+ONE ;
    real retval,PROD,DINDX,FRAC ;

    if (N <= 1) {
	PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA) ;
	PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA)) ;
	retval = PROD * pow(TWO,CONST)/(TWO*DN+CONST) ;
    } else {
	PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE) ;
	PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE)) ;
	PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA) ;
	PROD  = PROD*(ONE+BETA)*(TWO+BETA) ;
	for (int I=3 ; I<= N ; I++) {
	    DINDX = (real) I ;
	    FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA)) ;
	    PROD  = PROD*FRAC ;
	}
	retval = PROD * pow(TWO,CONST)/(TWO*DN+CONST) ;
    }
  
    return retval ;
}


void LAGF(RealVector* P,RealVector* PDER, RealVector* PDER2, RealVector* ZQ, RealVector* ZP, int M)
{
    // Computes the Lagrange polynomial P of degree 
    // NP-1 and its derivatives at ZQ 

    int  I, J, L, N ;
    int  IQ ;
    RealVector *PROD, *SUM ;

    // Initialize
    int  NP = ZP -> GetSize() ;
    int  NQ = ZQ -> GetSize() ;
    real COEFF = ONE ;

    PROD = new RealVector(NQ) ; SUM = new RealVector(NQ) ;

    for (I = 1 ; I <= NQ ; I++) {
	P     -> At(I) = ONE ;
	PDER  -> At(I) = ONE ;
	PDER2 -> At(I) = ONE ;
    }

    // Lagrange polynomial
    for (J = 1 ; J <= M-1 ; J++) {
	for (IQ = 1 ; IQ <= NQ ; IQ++) {
	    P->At(IQ) = P->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
	}
	COEFF = COEFF/(ZP->At(M)-ZP->At(J)) ;
    }
    for (J = M+1 ; J <= NP ; J++) {
	for (IQ = 1 ; IQ <= NQ ; IQ++) {
	    P->At(IQ) = P->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
	}
	COEFF = COEFF/(ZP->At(M)-ZP->At(J)) ;
    }
    P -> Multiply(COEFF) ;

    // First derivative of the Lagrange polynomial
    for (L = 1 ; L <= M-1 ; L++) {
	for (I = 1 ; I <= NQ ; I++) {
	    PROD->At(I) = ONE ;
	}
	for (J = 1 ; J <= NP ; J++) {
	    if ((J != L) && (J != M)) {
		for (IQ = 1 ; IQ <= NQ ; IQ++) {
		    PROD->At(IQ) = PROD->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
		}
	    }
	}
	PDER -> Add(PROD) ;
    }
    for (L = M+1 ; L <= NP ; L++) {
	for (I = 1 ; I <= NQ ; I++) {
	    PROD->At(I) = ONE ;
	}
	for (J = 1 ; J<= NP ; J++) {
	    if ((J != L) && (J != M)) {
		for (IQ = 1 ; IQ <= NQ ; IQ++) {
		    PROD->At(IQ) = PROD->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
		}
	    }
	}
	PDER -> Add(PROD) ;
    }
    PDER -> Multiply(COEFF) ;

    // Second derivative of the Lagrange polynomial
    for (L = 1 ; L <= M-1 ; L++) {
	for (I = 1 ; I <= NQ ; I++) {
	    SUM->At(I) = ZERO ;
	}
	for (N = 1 ; N <= NP ; N++) {
	    if ((N != L) && (N != M)) {
		for (I = 1 ; I <= NQ ; I++) {
		    PROD->At(I) = ONE ;
		}
		for (J = 1 ; J <= NP ; J++) {
		    if((J !=  L) && (J != M) && (J != N)) {
			for (IQ = 1 ; IQ <= NQ ; IQ++) {
			    PROD->At(IQ) = PROD->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
			}
		    }
		}
		SUM -> Add(PROD) ;
	    }
	}
	PDER2 -> Add(SUM) ;
    }
    for (L = M+1 ; L <= NP ; L++) {
	for (I = 1 ; I <= NQ ; I++) {
	    SUM->At(I) = ZERO ;
	}
	for (N = 1 ; N <= NP ; N++) {
	    if ((N != L) && (N != M)) {
		for (I = 1 ; I <= NQ ; I++) {
		    PROD->At(I) = ONE ;
		}
		for (J = 1 ; J <= NP ; J++) {
		    if ((J != L) && (J != M) && (J != N)) {
			for (IQ = 1 ; IQ <= NQ ; IQ++) {
			    PROD->At(IQ) = PROD->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
			}
		    }
		}
		SUM -> Add(PROD) ;
	    }
	}
	PDER2 -> Add(SUM) ;
    }
    PDER2 -> Multiply(COEFF) ;
}

void LAGF(RealVector* P, RealVector* ZQ, RealVector* ZP, int M)
{
    // Computes the Lagrange polynomial P of degree 
    // NP-1

    int  I, J, L, N ;
    int  IQ ;
    RealVector *PROD, *SUM ;

    // Initialize
    int  NP = ZP -> GetSize() ;
    int  NQ = ZQ -> GetSize() ;
    real COEFF = ONE ;


    for (I = 1 ; I <= NQ ; I++) {
	P->At(I) = ONE ;
    }

    for (J = 1 ; J <= M-1 ; J++) {
	for (IQ = 1 ; IQ <= NQ ; IQ++) {
	    P->At(IQ) = P->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
	}
	COEFF = COEFF/(ZP->At(M)-ZP->At(J)) ;
    }
    for (J = M+1 ; J <= NP ; J++) {
	for (IQ = 1 ; IQ <= NQ ; IQ++) {
	    P->At(IQ) = P->At(IQ)*(ZQ->At(IQ)-ZP->At(J)) ;
	}
	COEFF = COEFF/(ZP->At(M)-ZP->At(J)) ;
    }
    P -> Multiply(COEFF) ;
}
