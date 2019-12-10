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

#ifndef JACOBI
#define JACOBI
/*!
 * \file jacobi.hxx
 * \brief This file contains C++ versions of the Fortran routines of Einar Ronquist's MIT library
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/param.hxx"

class FullMatrix ;
class RealVector ;


void DGJD    (FullMatrix* d, RealVector* z, real alpha, real beta) ;
void DGLJD   (FullMatrix* d, RealVector* z, real alpha, real beta) ;

void IGJMD   (FullMatrix* i12, RealVector* z1, RealVector* z2, real alpha, 
	      real beta) ;
void IGLJMD  (FullMatrix* i12, RealVector* z1, RealVector* z2, real alpha, 
	      real beta, int interior=0) ;            

void ZIGLJD  (RealVector* z, real alpha, real beta) ;
void ZWGJD   (RealVector* z, RealVector* w, real alpha, real beta) ;
void ZGJD    (RealVector* z, real alpha, real beta) ;
void ZWGLJD  (RealVector* z, RealVector* w, real alpha, real beta) ;
void ZGLJD   (RealVector* z, real alpha, real beta) ;
void ZWIGLJD (RealVector* z, RealVector* zw, RealVector* w, real alpha, real beta) ;

real ENDW1   (int n, real alpha, real beta) ;
real ENDW2   (int n, real alpha, real beta) ;
real GAMMAF  (real x) ;
void JACG    (RealVector* xjac, real alpha, real beta) ;
void JACOBF  (real& poly, real& pder, real& polym1, real& pderm1, real& polym2,
	      real& pderm2, int n, real alp, real bet, real x) ;
real PNORMJ  (int n, real alpha, real beta) ;

void PHIML   (FullMatrix* mat, RealVector* z);  
real PNLEG   (real z, int n);

void LAGF    (RealVector* x1, RealVector* x2, RealVector* x3, RealVector* z1, RealVector* z2, int M) ;
void LAGF    (RealVector* x1, RealVector* z1, RealVector* z2, int M) ;

#endif
