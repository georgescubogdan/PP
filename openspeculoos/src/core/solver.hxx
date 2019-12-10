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

#ifndef SOLVER
#define SOLVER

/*!
 * \file solver.hxx
 * \brief Manages solvers
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/garbaged.hxx"
#include "core/vector.hxx"
#include "core/param.hxx"
#include "core/matrix.hxx"
#include "core/timer.hxx"
#include <string>

class SemiDiscreteProblem ;
class Field ;
// constant used to drive the solvers

const real     TINY      = 1.e-19 ;
const real     VERYSMALL = 1.e-15 ;
const real     SMALL     = 1.e-12 ;
const real     BIG       = 1.e+5 ;
const real     LARGE     = 1.e+7 ;
const real     VLARGE    = 1.e+10 ;

/*! \class Solver
   * \brief Base class managing solvers, both iterative and direct.
   * 
   * A solver is defined by its verbosity and the unknown field to
   * find. The problem to solve will be passed in argument to the
   * virtual method Solve.
   * 
   * \see {DirectSolver , IterativeSolver , Problem}
   */ 
class Solver : public GarbagedObject
{

protected :

    string name ;
    int    verbosity ;

public :

    // = CONSTRUCTORS

    Solver () ;
    virtual ~Solver () ;

    // = INQUIRY AND ATTRIBUTE SETTING

    virtual void SetAdjustTolerance (boolean b) ;
    void         SetName (string s) ;
    void         SetVerbosity (int v) ;

    // = SOLVE

    virtual void Solve (SemiDiscreteProblem* prob) = 0 ;

    // = PRINTINGS

    virtual void Print () ;
    virtual void PrintTime () ;
    virtual void PrintTime (FILE *fp) ;
    // debug

} ;

/*! \class ConjugateGradient
   * \brief Manages the preconditioned conjugate gradient algorithm.
   * 
   * The maximum number of iterations and the required tolerance are the
   * additional attributes.
    //
   * The algorithm is described in: Gene H. Golub and C.F. Van Loan, 
   * "Matrix Computations", 2nd edition, The Johns Hopkins University, 
   * p 374, 1983.
   * 
   * The convergence criterion is: norm(residue)/norm(RHS) smaller than tolerance,
   * where 'residue' and 'RHS' are equal at iteration 0.
   * 
   * The user may want to relax the convergence tolerance. This situation 
   * occurs if two conjugate-gradient solvers are nested: the user may
   * want the convergence tolerance of the inner solver to increase in the
   * same proportion as the outer solver converges, in order to avoid
   * requesting from the inner solver to keep its tolerance too demanding.
   * If <{adjustTolerance}> is 'false' (the default), then <{normRHS}> is
   * evaluated at the beginning of the iterative process, otherwise the
   * value of <{normRHS}> calculated in a previous iterative process is
   * used.
   * Every time the solver receives the message SetAjustTolerance(true),
   * the solver forgets the value of <{normRHS}>.
   * 
   * The user may want to request from the solver the number of performed
   * iterations for its last run (attribute <{iter}>), as well as the
   * cumulated number of iterations for all runs (attribute 
   * <{iterCumulated}>). The latter is useful when the solver is used
   * repeatedly (e.g., time stepping or nested solvers). Use method 
   * Print() for printing these numbers.
   */ 
class ConjugateGradient : public Solver 
{

protected :

    int     maxIterations ;
    real    tolerance ;
    real    normRHS ;
    boolean adjustTolerance ;
    int     iter ;
    int     iterCumulated ;

    // TIMER DEBUG EP
    float*   TimePCGtot ;
    float*   TimePCG_Ax ;

public :

    // = CONSTRUCTORS

    ConjugateGradient (int maxIter, real tol) ;
    ~ConjugateGradient () ;

    // = SETTING ;

    void SetAdjustTolerance (boolean b=false) ;

    // = SOLVE
 
    void Solve (SemiDiscreteProblem* prob) ;

    // = PRINTINGS

    void Print () ;
} ;

/*! \class GMRES
   * \brief Manages the preconditioned conjugate gradient algorithm.
   * 
   * The maximum number of iterations and the required tolerance are the
   * additional attributes.
   * 
   * The algorithm is described in: Gene H. Golub and C.F. Van Loan, 
   * "Matrix Computations", 2nd edition, The Johns Hopkins University, 
   * p 374, 1983.
   * 
   * The convergence criterion is: norm(RHS)/norm(residue) smaller than tolerance,
   * where 'residue' and 'RHS' are equal at iteration 0.
   * 
   * The user may want to relax the convergence tolerance. This situation 
   * occurs if two conjugate-gradient solvers are nested: the user may
   * want the convergence tolerance of the inner solver to increase in the
   * same proportion as the outer solver converges, in order to avoid
   * requesting from the inner solver to keep its tolerance too demanding.
   * If <{adjustTolerance}> is 'false' (the default), then <{normRHS}> is
   * evaluated at the beginning of the iterative process, otherwise the
   * value of <{normRHS}> calculated in a previous iterative process is
   * used.
   * Every time the solver receives the message SetAjustTolerance(true),
   * the solver forgets the value of <{normRHS}>.
   * 
   * The user may want to request from the solver the number of performed
   * iterations for its last run (attribute <{iter}>), as well as the
   * cumulated number of iterations for all runs (attribute 
   * <{iterCumulated}>). The latter is useful when the solver is used
   * repeatedly (e.g., time stepping or nested solvers). Use method 
   * Print() for printing these numbers.
   */ 
class GMRES : public Solver 
{

protected :

    int     maxIterations ;
    real    tolerance ;
    real    normRHS ;
    boolean adjustTolerance ;
    int     iter ;
    int     iterCumulated ;

    // TIMER DEBUG EP
    float*   TimePCGtot ;
    float*   TimePCG_Ax ;

public :

    // = CONSTRUCTORS

    GMRES (int maxIter, real tol) ;
    ~GMRES () ;

    // = SETTING ;

    void SetAdjustTolerance (boolean b=false) ;

    // = SOLVE
 
    void Solve (SemiDiscreteProblem* prob) ;

    // = PRINTINGS

    void Print () ;
} ;

/*! \class DirectSolver
   * \brief Manages direct solver based on Lapack
   * 
   * The matrix of the system is built by probing. Consequently, the
   * approach is very general but expensive. This class is thus a debugging
   * tool; the standard solver is the conjugate gradient one.
   * 
   * Optionally the user may request that the value of one unknown (chosen
   * as the last unknown) of the system be prescribed to 0. In that case
   * attribute <{suppress}> is set to 'true' (default: 'false'). This is
   * useful when pressure operators, which are singular, are concerned. 
   */ 

class DirectSolver : public Solver
{

protected :

    FullMatrix*   matrix;
    boolean       storeMatrix;
    boolean       suppress ;

public :

    // = CONSTRUCTORS

    DirectSolver (boolean store = false, boolean opt = false) ;

    // = SOLVE

    void Solve (SemiDiscreteProblem* prob) ;

} ;


#endif
