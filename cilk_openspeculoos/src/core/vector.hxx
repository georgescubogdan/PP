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

#ifndef VECTOR
#define VECTOR

/*!
 * \file vector.hxx
 * \brief Vector utilities
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/util.hxx"
#include "core/options.hxx"
#include "core/vector_templates.hxx"
#include <stdio.h>

/*! \class RealVector
   * \brief Manages vectors storing real numbers. 
   * 
   * BLAS routines are used whenever possible for speeding up calculations.
   */ 


class RealVector : public Vector<real>
{

public :

    // = CONSTRUCTORS

    inline RealVector () ; 
    inline RealVector (int n) ; 

    // = OPERATIONS

    inline void    Add (RealVector* x) ;
    inline void    Add (RealVector* x, real alpha) ;
    inline void    CopyFrom (RealVector* x) ;
    void           Divide (RealVector* x) ;
    real           Dot (RealVector* x) ;
    RealVector*    Duplicate () ;
    real           GetEuclidianNorm () ;
    real           GetInfiniteNorm () ;
    void           Inverse () ;
    void           Print () ;
    void           Add (real alpha) ;

    void           Multiply (real alpha) ;
    void           Multiply (RealVector* x) ;
    void           MultiplyAndAdd (real alpha, RealVector* x, real beta=ONE) ;
    void           MultiplyAndSwitchSigns (RealVector* x) ;
    void           SetValuesZero () ;
    void           SetToConstant (real) ;
    void           Square() ;
    void           Sqrt () ;
    void           Power(real x) ; // MAH 04.2006
    real           SumValues () ;
} ;

/*! \class Mpi_Buffer
   * \brief Manages temporary repositories for exchanging values between processors
   * 
   * Performing an MPI transaction every time a spectral element must send
   * its values to a spectral element assigned to another processor is 
   * unacceptably time consuming. An MPI buffer is an array that allows all
   * elements on a processor p to store the values they have to send to (or
   * to receive from) elements on another processor q. 
   * Thanks to MPI buffers, for any processor pair (p,q), instead of many
   * MPI transactions, only one MPI transaction is needed: for sending/
   * receiving the buffer (<{size}> real numbers).
   * 
   * An MPI buffer is used either for storing values (to be sent using MPI)
   * or for retrieving values (already received by using MPI):
   * . in the storing case, the elements append their values to the buffer
   *   using method Append, making the buffer size grow from 0 to its full
   *   size <{size}>. Then these <{size}> values will be sent to another
   *   processor q during the next processor communication phase;
   * . in the retrieving case, <{size}> values have been received from
   *   another processor q during the last processor communication phase.
   *   The elements read their values from the buffer, in the same order as
   *   they have been stored by the other processor, using method Get, and
   *   making the buffer index <{start}> grow from 0 to size-1.
   * 
   * Mpi buffers are not used on a scalar implementation.
   */ 


class Mpi_Buffer : public RealVector
{

protected :

    int start ;

public :

    // = CONSTRUCTORS

    Mpi_Buffer () ; 

    // = OPERATIONS

    void Append (int n, real* x, int incr) ;
    void Get (int n, real* x, int incr, int mode=0) ;
    int  GetStart () ;
    void Reset () ;
} ;


//___________________________ RealVector (inline) _____________________________


RealVector  :: RealVector ()
    : Vector<real> ()
{
    // Constructor. Initializes the receiver as a vector of size 0.
}


RealVector  :: RealVector (int n)
    : Vector<real> (n)
{
    // Constructor. Initializes the receiver as a vector of size 'n' with all
    // values initialized to 0.

    SetValues(ZERO) ;
}


void RealVector :: Add (RealVector* x)
{
    // y = y + x .

    Add(x,ONE) ;
}


void RealVector :: Add (RealVector* x, real alpha)
{
    // y = y + alpha x .
    // Adds the product 'alpha*x' to the receiver 'y'.

# ifdef REQUIRE
    Require("y and x have same size", size == x->size) ;
    InMethod("RealVector::Add(x,alpha)") ;
# endif

    int ione=IONE ;

    DAXPY(&size, &alpha, &x->values[0], &ione, &values[0], &ione) ;
}


void RealVector :: CopyFrom (RealVector* x)
{
    // y = x .
    // Copies 'x' into the receiver 'y'.

# ifdef REQUIRE
    Require("y and x have same size", size == x->size) ;
    InMethod("RealVector::CopyFrom(x)") ;
# endif

    int ione=IONE ;

    DCOPY(&size, &x->values[0], &ione, &values[0], &ione) ;
}


#endif

