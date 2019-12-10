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

#ifndef UTIL
#define UTIL

/*!
 * \file util.hxx
 * \brief Utilities
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/param.hxx"
#include <string>


class Point ;


// This file implements utility functions.


// --------------
// ERROR MESSAGES
// --------------

void  Error (string method, string message) ;
void  ImplementedBySubclasses (string method) ;
void  InMethod (string method) ;
void  NotImplemented (string method) ;
void  NotImplemented (string method, int flag) ;
void  Require   (string message, boolean test) ;
void  Warning (string method, string message) ;
void  Exit (int code=0) ;


// ------------------
// DYNAMIC ALLOCATION
// ------------------

void  ExhaustedMemory () ;


// -------------
// INTERPOLATION
// -------------

real  Phi1 (real r) ;
real  Phi2 (real r) ;
real  ParentToPhysic (real x1, real x2, real r) ;


// ----------------------
// MATHEMATICAL FUNCTIONS
// ----------------------

real  EvalLegendre (int n, real x) ;
int   Factorial (int n) ;
real  Pow (real x, int n) ;

real  angl2pi(real COS,real SIN) ;

// -------------
// MISCELLANEOUS
// -------------

int   CheckConvexity (Point*, Point*, Point*, Point*) ;
FILE* Fopen (const char* fileName, const char* mode) ;

// -----------------------
// ONE-PROCESSOR PRINTINGS
// -----------------------

#define Printf if(Mpi_rank==0) printf
#define FPrintf if(Mpi_rank==0) fprintf
// use of this macro requires including file parallel.hxx


#endif
