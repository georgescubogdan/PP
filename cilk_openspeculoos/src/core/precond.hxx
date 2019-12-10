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

#ifndef PRECOND
#define PRECOND

/*!
 * \file precond.hxx
 * \brief Manages preconditionners
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/vector.hxx"
#include "core/garbaged.hxx"

class SemiDiscreteProblem ;
class Field ;
class TensorMatrix ;
class FullMatrixLU ;
class DiagonalMatrix ;

enum CouzyType { DIRICHLET_NEUMANN, DIRICHLET } ;


/*! \class Preconditioner
   * \brief Manages post-processors output.
   * 
   * A preconditioner P acts on a <{problem}> which implements a linear 
   * system A x = rhs.
   * The task of the preconditioner is to solve P z = r, with P supposedly
   * "close" to A.
   */ 

class Preconditioner : public GarbagedObject
{ 

protected :

    SemiDiscreteProblem* problem ;

public :

    // = CONSTRUCTORS

    Preconditioner () ;
    virtual ~Preconditioner () ;

    // = SETTINGS
 
    virtual void SetProblem (SemiDiscreteProblem* prob) ;

    // = PRECONDITIONING

    virtual void Precondition (Field* r, Field* z) = 0 ;
} ;


/*! \class FieldPreconditioner
   * \brief Base class of preconditioners which can be stored in a simple field.
   * 
   * Solves P z = r.
   * The inverse of P is stored in <{inverseP}>.
   */ 

class FieldPreconditioner : public Preconditioner
{

protected :

    Field* inverseP ;
 
    virtual Field* GetInverseP () = 0 ;

public :

    // = CONSTRUCTORS

    FieldPreconditioner () ;
    ~FieldPreconditioner () ;

    // = SETTINGS

    void SetProblem (SemiDiscreteProblem* prob) ;

    // = PRECONDITIONING

    virtual void Precondition (Field* r, Field* z) ;
} ;


/*! \class MassPreconditioner
   * \brief Base class of preconditioners which can be stored in a simple field.
   * 
   * Solves P z = r, where P = B, the assembled mass matrix of the problem.
   */ 

class MassPreconditioner : public FieldPreconditioner
{

    friend class DirectSolver;
protected :

    Field* GetInverseP () ;

public :

    // = PRECONDITIONING

    void Precondition (Field* r, Field* z) ;
} ;


/*! \class DiagonalPreconditioner
   * \brief Diagonal preconditioner.
   * 
   * Solves P z = r, where P = diag(A).
   * The inverse of P, with the values equal to 1./0. replaced by 1., is 
   * stored in <{inverseP}>.
   */ 

class DiagonalPreconditioner : public FieldPreconditioner
{

protected :

    Field* GetInverseP () ;
} ;


/*! \class EuclidianNormPreconditioner
   * \brief Preconditioner by the euclidian norm of every column.
   * 
   * Solves P z = r, where P is the diagonal matrix such that, for any i, 
   * P(i,i) is the euclidian norm of the i-th column of A.
   */ 
class EuclidianNormPreconditioner : public FieldPreconditioner
{

protected :

    Field* GetInverseP () ;
} ;


/*! \class InfiniteNormPreconditioner
   * \brief Preconditioner by the infinite norm of every column.
   * 
   * Solves P z = r, where P is the diagonal matrix such that, for any i, 
   * P(i,i) is the infinite norm of the i-th column of A.
   */ 

class InfiniteNormPreconditioner : public FieldPreconditioner
{

protected :

    Field* GetInverseP () ;
} ;


/*! \class PseudoLaplacePreconditioner
   * \brief Element-by-element preconditioner for a pseudo-Laplace operator.
   * 
   * Solves P z = r, where P is the pseudo-Laplace operator
   * P = E = D inv(B) D(T)
   * on a pressure (i.e., GL) grid.
   * 
   * E is inverted approximately, element by element.
   * 
   * A particular problem is the kind of boundary conditions to consider on
   * every element for the inversion of block(B) in block(E). A part C_e of
   * the contour of an element is either:
   * - on the contour of the mesh subject to Dirichlet conditions. In
   * this case, homogenous Dirichlet conditions are applied to block(B)
   * on C_e;
   * - on the contour of the mesh subject to Neumann conditions. In this
   * case, homogenous Neumann conditions (i.e., nothing) are applied to
   * block(B) on C_e;
   * - on the interface with another element. In this case,
   * . if attribute <{type}> is set to DIRICHLET_NEUMANN, homogeneous
   * Neumann conditions (i.e., nothing) are applied to block(B) on
   * C_e. This is the default.
   * . if attribute <{type}> is set DIRICHLET, homogenous Dirichlet
   *    conditions are applied to block(B) on C_e. 
   * 
   * W. Couzy advocates for the DIRICHLET_CONDITION case, which is the
   * default. The DIRICHLET 
   * 
   * On block(E) no boundary conditions are considered, since E is defined
   * on a pressure (GL) grid.
   *  
   * Attribute <{elementaryMatrices}> stores matrix E of every element,
   * including the boundary conditions. The elementary systems E z = r are
   * solved by LU decomposition using Lapack.
   * 
   * References: see CouzyPreconditioner.
   */ 

class PseudoLaplacePreconditioner : public Preconditioner
{

protected :

    CouzyType              type ;
    Vector<FullMatrixLU*>* elementaryMatrices ;

    void SetElementaryMatrices () ;

public :

    // = CONSTRUCTORS

    PseudoLaplacePreconditioner (CouzyType=DIRICHLET_NEUMANN) ;
    ~PseudoLaplacePreconditioner () ;

    // = PRECONDITIONING

    void Precondition (Field* r, Field* z) ;
} ;


/*! \class CouzyPreconditioner
   * \brief Couzy's preconditioner for (Navier-)Stokes problems.
   * 
   * Solves P z = r, where P is the P4 preconditioner deviced by Wouter
   * Couzy.
   * 
   * This preconditioner can be used for a problem of type SteadyStokes, in
   * either Uzawa, BP1 or BP1-PC decomposition mode.
   * 
   * P is an approximation to "minus" the Uzawa or BP1 operator 
   *     S = D inv(H) D(T),
   * where:
   *   . H = lambda B - nu L (case Uzawa);
   *   . H = lambda B        (case BP1 or BP1-PC);
   *   . D is the weak divergence operator;
   *   . B is the weak identity (mass) matrix on the velocity grid;
   *   . L is the (negative definite) weak Laplacian operator;
   *   . lambda and nu are nonnegative constants.
   * It relies on the assumption that the inverse of S can be approximated
   * (modulo a multiplicative constant) by the inverse of P defined as:
   *     inv(P) = nu inv(B~) + lambda inv(E)  (case Uzawa), 
   *     inv(P) =                     inv(E)  (case BP1 or BP1-PC), 
   * where:
   *   . B~ is the mass matrix on the pressure grid. B~ is diagonal and so
   *     is trivially inverted;
   *   . E = D inv(B) D(t), a pseudo-Laplace operator. Since E is as 
   *     difficult to invert as S, block(E) is considered instead of E. 
   *     block(E) is the element-by-element (i.e., unassembled)
   *     approximation to E. block(E) is inverted by a direct method (LU)
   *     on every element.
   * 
   * Attribute <{pseudoLaplace}> performs the inversion of block(E). See in
   * class PseudoLaplacePreconditioner the matter related to boundary
   * conditions.
   * 
   * The derivation of this preconditioner can be found in either of the
   * following equivalent references:
   *   [1] W. Couzy and M.O. Deville, Spectral-Element Preconditioners
   *       for the Uzawa Pressure Operator Applied to Incompressible Flows,
   *       J. of Scientific Computing, Vol. 9, No. 2, pp 107-122, 1994.
   *   [2] W. Couzy, Spectral Element Discretization of the Unsteady
   *       Navier-Stokes Equations and its Iterative Solution on Parallel
   *       Computers, Ph.D. thesis No. 1380, DGM-LMF, Swiss Federal
   *       Institute of Technology, 1015 Lausanne, Switerland, 1995.
   */ 

class CouzyPreconditioner : public Preconditioner
{
 
protected :

    Preconditioner* pseudoLaplace ;

public :

    // = CONSTRUCTORS

    CouzyPreconditioner (CouzyType=DIRICHLET_NEUMANN) ;
    ~CouzyPreconditioner () ;

    // = SETTINGS
 
    void SetProblem (SemiDiscreteProblem* prob) ;

    // = PRECONDITIONING

    void Precondition (Field* r, Field* z) ;
} ;


#endif
