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

#ifndef PROBLEM
#define PROBLEM

/*!
 * \file problem.hxx
 * \brief Manages problem types
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/vector.hxx"
#include "core/garbaged.hxx"
#include "core/param.hxx"

#define cout_root if(Mpi_rank==0) std::cout

class SemiDiscreteProblem ;
class Helmholtz ;
class SteadyStokes ;
class Mesh ;
class DirichletCondition ;
class NeumannCondition ;
class Solver ;
class Preconditioner ;
class TimeIntegrationScheme ;
class Field ;
class FlatField ;
class Function ;
class MortaredField ;
class FullMatrix ;

enum DecompositionType { UZAWA, BP1, BP3, BP1_PC, BP3_PC, YOSIDA } ;

/*! \class Problem
   * \brief Base class for all types of problems.
   * 
   * A problem is a differential equation on a domain, submitted to 
   * Dirichlet and Neumann boundary conditions.
   * The domain is implemented by <{mesh}>.
   * For economy, a forcing-term field <{forcing}> and a function for the
   * forcing-term field are implemented at this level.
   * 
   * The functional F often contains a forcing-term f(t); moreover, the 
   * prescribed values of the Dirichlet conditions may be time dependent.
   * Therefore the current value of t is stored in <{time}>.
   */ 
class Problem : public GarbagedObject
{ 

protected :

    Mesh*                        mesh ;
    Vector<DirichletCondition*>* dirichletConditions ;       // -> SDP
    Vector<NeumannCondition*>*   neumannConditions ;         // -> SDP
    FlatField*                   forcing ;                   //   id.
    Function*                    forcingFunction ;           //   id.
    real                         time ;                      //   id.

    void DeleteDirichletConditions () ;
    void DeleteNeumannConditions () ;

	static FILE* ficdeb; 

public :

    // = CONSTRUCTORS

    Problem (Mesh* m) ;
    virtual ~Problem () ;

    // = INQUIRY

    Mesh*                     GetMesh () ;

    // = DIRICHLET CONDITIONS

    virtual void              AddDirichletCondition (DirichletCondition* dc) ;
    virtual void              SetDirichletConditions 
    (Vector<DirichletCondition*>* dcs) ;
    void                      ApplyDirichletConditions (Field* x) ;
    void                      ApplyHomogeneousDirichletConditions (Field* x) ;
    void                      ApplyHomogeneousDirichletConditions (Field* x,
                                                                   int icx) ;
    Vector<DirichletCondition*>* GetDirichletConditions () ;
    Vector<DirichletCondition*>* GetHomogeneousDirichletConditions () ;

    // = NEUMANN CONDITIONS
    virtual Vector<NeumannCondition*>* GetNeumannConditions();
    virtual void              AddNeumannCondition (NeumannCondition* nc) ;
    virtual void              SetNeumannConditions 
    (Vector<NeumannCondition*>* ncs) ;

    // = FORCING TERM
 
    virtual void              SetForcing (FlatField* f) ;
    virtual void              SetForcingFunction (Function* funct) ;
    void                      EvaluateForcing (real t) ;
    FlatField*                GetForcing () ;

    // = TIME

    virtual void              SetTime (real t) ;
} ;


/*! \class HeatEquation
   * \brief Solves a heat-equation problem.
   * 
   * The problem consists in the parabolic partial-differential equation
   * 
   * du/dt = F(u,t)   (+ Dirichlet, Neumann and initial conditions)
   *                 2
   * with F(u,t) = C  Laplacian(u) + f(x,t).
   *      2                                           2
   * - C can be viewed as the thermal diffusivity: C = k / (rho c), where
   * k is the Fourier coefficient, rho the mass density and c the 
   * specific heat;
   * - the forcing term f can be viewed as a source term (power per unit
   * volume) divided by rho*c.
   * 
   * The problem is solved step by step, by means of its attribute 
   * <{timeIntegrationScheme}>. At every step this scheme solves a time-
   * discretized problem of type Helmholtz.
   */ 

class HeatEquation : public Problem
{

protected :

    Field*                 u ;
    TimeIntegrationScheme* timeIntegrationScheme ;
    Helmholtz*             helmholtz ;

public :

    // = CONTRUCTORS

    HeatEquation (Mesh* m) ;
    ~HeatEquation () ;

    // = SETTINGS AND ACCESSINGS

    void       AddDirichletCondition (DirichletCondition* dc) ;
    void       AddNeumannCondition (NeumannCondition* nc) ;
    Helmholtz* GetHelmholtz () ;
    void       SetField (Field* x) ;
    void       SetForcing (FlatField* f) ;
    void       SetForcingFunction (Function* funct) ;
    void       SetIntegrationRule (FlatField* rule) ;
    void       SetPreconditioner (Preconditioner* p) ;
    void       SetSolver (Solver* s) ;
    void       SetThermalDiffusivity (real val) ;
    void       SetTime (real t) ;
    void       SetTimeIntegrationScheme (TimeIntegrationScheme* tis) ;

    // = SOLUTION

    void       Solve (int nstep) ;
    Field*     GetUnknown () ;
} ;


/*! \class Stokes
   * \brief Solves an unsteady incompressible Stokes problem.
   * 
   * The problem consists in the coupled PDE system
   * 
   * dv/dt = - grad p + nu Laplacian v + f - div v = 0
   * (+ Dirichlet, Neumann and initial conditions on v).
   * 
   * The problem is solved step by step, by means of its attribute 
   * <{timeIntegrationScheme}>. At every step this scheme solves a time-
   * discretized problem of type SteadyStokes.
   */ 
class Stokes : public Problem
{

protected :

    MortaredField*         velocity ;
    FlatField*             pressure ;
    TimeIntegrationScheme* timeIntegrationScheme ;
    SteadyStokes*          steadyStokes ;

public :

    // = CONTRUCTORS

    Stokes (Mesh* m, DecompositionType=UZAWA) ;
    virtual ~Stokes () ;

    // = ATTRIBUTE SETTING AND ACCESSING

    void           AddDirichletConditionV (DirichletCondition* dc) ;
    void           AddNeumannConditionV (NeumannCondition* nc) ;
    SteadyStokes*  GetSteadyStokes () ;
    void           SetForcing (FlatField* f) ;
    void           SetForcingFunction (Function* funct) ;
    void           SetIntegrationRuleMass (FlatField* rule) ;
    void           SetIntegrationRuleMomentum (FlatField* rule) ;
    void           SetNu (real val) ;
    void           SetTimeIntegrationScheme (TimeIntegrationScheme* tis) ;
    virtual void   SetVelocity (MortaredField* v) ;
    void           SetPressure (FlatField* p) ;
    void           SetSolver (Solver* s) ;
    void           SetPreconditioner (Preconditioner* p) ;

    void           ImposeNoPressureNullSpace();
    // = SOLUTION

    void           Solve (int nstep) ;
    void           SetValues (int i, FlatField* x);

    MortaredField* GetVelocity () ;
    FlatField*     GetPressure () ;
} ;


/*! \class NavierStokes
   * \brief Solves an unsteady incompressible Navier-Stokes problem.
   * 
   * The problem consists in the coupled PDE system
   * 
   * dv/dt + v.grad v = - grad p + 1/Re Laplacian v + f
   *  \- div v            = 0
   * (+ Dirichlet, Neumann and initial conditions).
   * 
   * The treatment of the nonlinear convective term (v.grad v) is splitted
   * from the treatment of the linear term (Stokes part). An explicit
   * scheme is used for the former, wheareas an implicit or explicit scheme
   * is used for the latter.
   * 
   * The inherited attribute <{timeIntegrationScheme}> deals with the
   * Stokes part of problem. Attribute <{convectionScheme}> deals with the
   * convective part; at present it must be an explicit scheme.
   * 
   * Attribute <{nonlin}> stores the current value of -v.grad v.
   */ 
class NavierStokes : public Stokes
{

protected :

    TimeIntegrationScheme* convectionScheme ;
    FlatField*             nonlin ;

    void GetNonlinearTerm (FlatField* y) ;

public :

    // = CONSTRUCTORS

    NavierStokes (Mesh* m, DecompositionType=UZAWA) ;
    ~NavierStokes () ;

    // = SETTINGS

    void                   SetConvectionScheme (TimeIntegrationScheme* tis) ;
    void                   SetReynolds (real re) ;
    void                   SetVelocity (MortaredField* v) ;
    void                   SetValues (int i, FlatField* x) ;

    // = SOLUTION

    //    void                   Solve (int nstep) ;
    real                  Solve (int nstep) ;
} ;

/*! \class SemiDiscreteProblem
   * \brief Base class for time-independent problems.
   * 
   * Discretizes and solves for u(x) the weak form of:
   * lambda * u = F(u) + Hist   (plus Dirichlet and Neumann conditions),
   * with:
   *    . u = u(x) = the unknown field;
   *    . F(u)     = a differential operator on u, implemented in the 
   *                 subclasses;
   *    . Hist     = an optional "history" vector;
   *    . lambda   = a scalar.
   * 
   * The corresponding attributes are 
   *    . <{lambda}>;
   *    . <{unknown}>, for u, which can be either a flat field or a 
   *      mortared field;
   *    . <{integrationRule}>, to perform numerical quadratures;
   *    . <{strongistory}>, to store the part of Hist that is given in
   *      strong from and thus will have to be put into weak form;
   *    . <{weakHistory}>, to store the part of Hist that is given in weak
   *      form;
   * 
   * A semi-discrete problem is able to evaluate:
   *    . its operator F(u), in a weak form;
   *    . its operator A(u) = lambda u - LHS(u), where LHS(u) is the part
   *      of F(u) that depends on u (also in a weak form);
   *    . its right-hand side RHS, the part of F(u) that does not depend on
   *      u (also in a weak form).
   * 
   * If 'unknown' is interpolated with rules that do not extend until the
   * borders (typically, Gauss-Legendre rules, as opposed to Gauss-Lobatto-
   * Legendre rules), then attribute 'dirichletConditions' must be left
   * empty. This is typically the case when 'unknown' is a pressure field.
   * 
   * A semi-discrete problem is often used as the result of the time-
   * integration of a PDE problem. In that case, depending on the time-
   * integration scheme used, the semi-discrete problem must be solved
   * either in an implicit mode (the default) or in an explicit, 
   * simplified, mode. Further details on this distinction are implemented
   * by subclasses. The explicit/implicit mode is governed by attribute
   * <{mode}>.
   */ 
class SemiDiscreteProblem : public Problem
{
 
protected :

    real            lambda ;
    Field*          unknown ;
    FlatField*      integrationRule ;
    Solver*         solver ;
    FlatField*      weakHistory ;
    FlatField*      strongHistory ;
    int             mode ;
    Preconditioner* preconditioner ;
    Field*          inverseB ;

    enum {IMPLICIT, EXPLICIT} ;

public :

    // = CONSTRUCTORS

    SemiDiscreteProblem (Mesh* m) ;
    virtual ~SemiDiscreteProblem () ;

    // = SETTINGS

    virtual void     SetWeakHistory (FlatField* x) ;
    virtual void     SetStrongHistory (FlatField* x) ;
    virtual void     SetLambda (real val) ;
    void             SetIntegrationRule (FlatField* rule) ;
    void             SetPreconditioner (Preconditioner* p) ;
    void             SetSolver (Solver* s) ;
    virtual void     SetTemporalUnknown (Field* u) ;
    virtual void     SetToExplicit () ;
    virtual void     SetToImplicit () ;
    virtual void     SetUnknown (Field* u) ;

    // = ACCESSING

    FlatField*       GetIntegrationRule () ;
    real             GetLambda () ;
    Preconditioner*  GetPreconditioner () ;
    Solver*          GetSolver () ;
    Field*           GetUnknown () ;
    // tempo debug
    inline FlatField*       GetWeakHistory () {return weakHistory;}
    // = INQUIRY

    virtual boolean  IsHelmholtz () ;
    virtual boolean  IsSteadyStokes () ;

    // = LHS OPERATOR

    virtual void     GetAx (Field* u, Field* Au) = 0 ;
    virtual void     GetFlatAx (FlatField* u, FlatField* Au) = 0 ;

    // = RHS OPERATOR

    virtual void     GetRHS (Field* rhs) = 0 ;

    // = SOLVE

    virtual void     Solve () ;

    // = OTHER OPERATORS

    virtual void     ApplyInverseB (Field* y) ;
    void             DeleteInverseB () ;
    virtual void     GetFunctional (FlatField* x, real t, FlatField* y) = 0 ;
    Field*           GetInverseDiagonal () ;
    Field*           GetInverseLumpedDiagonal (int iop) ;
    virtual Field*   GetInverseB () ;
    FullMatrix*      GetMatrixOperator () ;

    // PRECONDITIONING (COUZY)

    virtual real     GetC1 () ;
    virtual real     GetC2 () ;
} ;

/*! \class Helmholtz
   * \brief Weak-form solver for a Helmholtz problem.
   * 
   * Discretizes and solves for u(x) the weak form of: 
   * lambda * u = F(u) + Hist  (plus Dirichlet and Neumann conditions), 
   * with:
   *    . F(u)          = alpha * Laplace u + f;
   *    . Laplace u     = div(grad u);
   *    . f = f(x)      = forcing term; 
   *    . Hist          = a "history" vector;
   *    . lambda, alpha = scalars (usually both positive).
   * 
   * The problem is solved in its weak form:
   *    lambda (u, psi) + alpha (grad u, grad psi) = (f+g+Hist, psi)  (1)
   * where 'g' is the rhs of the Neumann conditions grad(u) = g.
   * 
   * After discretization (with, say, Hist=0), this yields:
   *    H u = B f + B_contour g  (plus Dirichlet Conditions), 
   * where:
   *    . H = lambda B - alpha L,
   *    . B is the discrete weak identity operator (= mass matrix);
   *    . B_contour is the mass matrix on the contour;
   *    . L is the (negative definite) discrete weak Laplacian operator;
   *    . u and f are the discrete values of u(x) and f(x).
   * 
   * The problem can be solved in any of the following two contexts:
   *    . as a true Helmholtz problem: Hist is not used;
   *    . as a Helmholtz problem obtained by the time-discretization of a
   *      PDE problem such as a heat problem or a Stokes problem. This
   *      justifies the presence of an additional history vector in the 
   *      right-hand side.
   * 
   * Moreover, if the semi-discrete problem is the by-product of the time
   * discretization of a PDE problem, two events can occur:
   *    . if the time-integration scheme is explicit (and then attribute
   *      <{equationWriting}> is set to EXPLICIT_EW), the left-hand side of
   *      (1) is replaced by: lambda (u,psi).
   *      Note that although this "explicit" form does not require any
   *      inversion of linear system, it is not trivial, since it features
   *      both mortar constraints and Dirichlet conditions;
   *      PAS A JOUR.
   * Default values: <{equationWriting}> = IMPLICIT_EW,
   * 
   * Note that if lambda=0, alpha=-1 and Hist=0, a Poisson problem is 
   * obtained.
   */ 
class Helmholtz : public SemiDiscreteProblem
{

protected :

    real alpha ; 
    boolean  kioPressure;

public :

    // = CONSTRUCTORS

    Helmholtz (Mesh* m, boolean kioP=false) ;
    ~Helmholtz () ;

    // = SETTINGS

    void    SetAlpha (real val) ;
    void    SubtractWeakHistory (FlatField* hist) ;

    // = ACCESSING

    real    GetAlpha () ;

    // = INQUIRY

    boolean IsHelmholtz () ;

    // = OPERATORS

    void    GetAx (Field* u, Field* Au) ;
    void    GetRHS (Field* rhs) ;
    void    GetFlatAx  (FlatField* u, FlatField* Au) ;
    void    GetFlatRHS (FlatField* rhs) ;
    void    GetFunctional (FlatField* x, real t, FlatField* y) ;

    // PRECONDITIONING (COUZY)

    real    GetC1 () ;
    real    GetC2 () ;

    // DEBUG
    void Print();
} ;


/*! \class SteadyStokes
   * \brief Discretization a steady-state incompressible Stokes problem in weak form
   * 
   * Discretization a steady-state incompressible Stokes problem in weak
   * form, using Blair Perot's generalized LU formulation.
   * 
   * Solves: lambda v - nu Laplace v + grad p = f
   *          - div v                         = 0
   *                       (plus Dirichlet and Neumann conditions on v)
   * 
   * in weak form :
   *    lambda (v,psi) + nu (grad v, grad psi) - (p, grad psi) = (f+g,psi)
   *    - (div v, phi)                                         = 0
   * 
   * For a truly steady Stokes problem, lambda = 0.
   * For a steady Stokes problem resulting from the time discretization of
   * an unsteady Stokes or Navier-Stokes problem, lambda incorporates the
   * time stepping and the density, whereas f incorporates the history
   * vector.
   * 
   * The problem is reformulated as a Generalized block LU decomposition
   * (see Wouter Couzy's thesis, chapter 5).
   * Three decomposition types are proposed: Uzawa (exact decomposition),
   * Blair Perot 1 (first-order approximation) and Blair Perot 3 (third-
   * order approximation). Only Uzawa can be used if lambda = 0.
   * 
   * The discretized system can be written as:
   *      H v - D(T) p = B f + B_contour g   (+ Dirichlet conditions)
   *     -D v          = 0   
   * where lambda B - nu L (B = weak identity operator, L = weak Laplacian
   * operator, L is negative definite) and (T) stands for a transposition. 
   * 
   * The generalized block LU decomposition approximates this system by the
   * following system:
   * 
   *      H v - H J D(T) p = B f + B_contour g  (+ Dir. cond.)     (1)
   *     -D v              = 0
   * 
   * where J is an arbitrary matrix that determines the projection method:
   * Decomposition type:
   *    - Uzawa:            J = H-1
   *    - Blair Perot 1:    J = (lambda B)-1   (i.e., 'fractional step')
   *    - Blair Perot 1 PC: J = idem
   *    - Blair Perot 3:    J = (5.12, p. 80 Couzy's thesis) 
   *    - Blair Perot 3 PC: J = idem
   * 
   * The system (1) is solved in three steps:
   * 
   *   1)  H v~       = B f + B_contour g  (+ nonhomogeneous Dc's)
   *   2) -D J D(T) p = D v~               (+ homogeneous Dc's on v)
   *   3)  v          = v~ + J D(T) p      (            id.        )
   * 
   * In the pressure-correction (PC) mode is used, the right-hand side B f
   * is modified into B f + D(T) p and the unknown of equation 2 is
   * modified into
   *      -D J D(T) (p-pN) = D v~      (+ homogeneous bc's on v)
   * where pN is the pressure at the previous time step.
   *    
   * The user must supply attributes <{unknown}> (the inherited pressure 
   * field p) and <{velocity}>. The problem initializes their values. 
   * Since interface continuity must be enforced for the velocity,
   * <{velocity}> must be a mortared field, whereas <{unknown}> (the
   * pressure) must be a flat field.
   * In the case of a PC decomposition mode, attribute <{previousPressure}>
   * stores pN; if the user does not supply it, the problem initializes
   * its values to 0.
   * 
   * The user may or may not supply a <{forcing}> term f. If supplied, that
   * field must have the same number of components as the velocity. It may 
   * have different interpolation rules from those of the velocity.
   * 
   * The inherited attribute <{integrationRule}> is a 1-component flat
   * field whose rules are used for numerical quadratures of the mass 
   * equations. It is typically set to a duplicate of the pressure field.
   * Attribute <{integrationRuleMomentum}> plays the same role for the
   * momentum-conservation equation. It is typically set to a 1-component
   * duplicate of the [main field of the] velocity field.
   * 
   * The inherited attributes <{dirichletConditions}> and
   * <{neumannConditions}> contain boundary conditions for the pressure
   * field. Since b.c.'s on the pressure field are not implemented at
   * present, these 2 attributes are always empty.
   * 
   * When a Dirichlet or Neumann condition on the velocity condition is
   * assigned to a steady Stokes problem, the problem hands it over to its
   * Helmholtz problem.
   * 
   * An attribute <{helmholtz}>, which acts on the velocity, is used in
   * phase 1 of the generalized-LU solving process. If the Uzawa 
   * decomposition is chosen, another subproblem <{internalHelmholtz}> is
   * used for solving systems H x = y in the internal loop of phase 2 and 
   * in phase 3.
   * 
   * Attribute <{vTilde}> is a debugging object which takes a record of v~
   * in phase 1. It is normally not used.
   */ 

class SteadyStokes : public SemiDiscreteProblem
{

protected :

    DecompositionType type ;
    real              nu ; 
    MortaredField*    velocity ;
    FlatField*        previousPressure ;
    FlatField*        integrationRuleMomentum;
    Helmholtz*        helmholtz ;
    Helmholtz*        internalHelmholtz ;
    FlatField*        vTilde ;

    void CreateHelmholtz () ;
    void CreateInternalHelmholtz () ;

    // flag in case of essential bc on pressure (--> no null space)
    // default is false
    boolean            stressfreeBC;  

public :

    // = CONSTRUCTORS

    SteadyStokes (Mesh* m, DecompositionType dt) ;
    ~SteadyStokes () ;

    // = PARAMETERS

    real                      GetNu () ;
    void                      SetLambda (real val) ;
    void                      SetNu (real val) ;

    // = INQUIRY

    boolean                   IsSteadyStokes () ;

    // = FORCING TERM

    void                      SetForcing (FlatField* f) ;
    void                      SetForcingFunction (Function* funct) ;

    // = PRESSURE FIELD

    inline FlatField*         GetPressure () ;
    void                      SetPressure (FlatField* p) ;

    // for re strating
    FlatField*                GetPreviousPressure();
    // = VELOCITY FIELD

    inline MortaredField* GetVelocity () ;
    void                  SetVelocity (MortaredField* v) ;
    void                  SetTemporalUnknown (Field* v) ;
    FlatField*            GetVtilde ()                     {return vTilde;}

    // = PRESSURE DIRICHLET CONDITIONS (NOT IMPLEMENTED)

    void                  AddDirichletCondition (DirichletCondition* dc) ;
    void                  SetDirichletConditions
    (Vector<DirichletCondition*>* dcs) ;

    // = PRESSURE NEUMANN CONDITIONS (NOT IMPLEMENTED)

    void                  ImposeNoPressureNullSpace();

    void                  AddNeumannCondition (NeumannCondition* nc) ;
    void                  SetNeumannConditions
    (Vector<NeumannCondition*>* ncs) ;

    // = VELOCITY DIRICHLET CONDITIONS

    void                  AddDirichletConditionV (DirichletCondition* dc) ;
    void                  ApplyHomogeneousDirichletConditionsV (Field *x) ;
    Vector<DirichletCondition*>* GetHomogeneousDirichletConditionsV () ;
    void                  SetDirichletConditionsV
    (Vector<DirichletCondition*>* dcs) ;

    // = VELOCITY NEUMANN CONDITIONS

    void                  AddNeumannConditionV (NeumannCondition* nc) ;
    void                  SetNeumannConditionsV
    (Vector<NeumannCondition*>* ncs) ;

    // = INTEGRATION RULES

    FlatField*            GetIntegrationRuleMass () ;
    FlatField*            GetIntegrationRuleMomentum () ;
    void                  SetIntegrationRuleMass (FlatField* rule) ;
    void                  SetIntegrationRuleMomentum (FlatField* rule) ;

    // = GENERALIZED LU DECOMPOSITION

    void                  ApplyJ (MortaredField* x, MortaredField* y) ;
    Helmholtz*            GetHelmholtz () ;
    Helmholtz*            GetInternalHelmholtz () ;

    // = OPERATORS

    void                  GetAx (Field* x, Field* y) ;
    void                  GetFlatAx (FlatField* x, FlatField* y) ;
    void                  GetRHS (Field* rhs) ;
    void                  GetFunctional (FlatField* x, real t, FlatField* y) ;
    void                  SetWeakHistory (FlatField* x) ;
    void                  SetStrongHistory (FlatField* x) ;

    FullMatrix*           GetMatrixOperatorV(int q, real coeff );
    // = IMPLICIT/EXPLICIT MODE

    void                  SetToExplicit () ;
    void                  SetToImplicit () ;

    // = TIME

    void                  SetTime (real t) ;

    // = SOLVE

    void                  Solve () ;

    // PRECONDITIONING (COUZY)

    real                  GetC1 () ;
    real                  GetC2 () ;
} ;


//_____________________________ SteadyStokes (inline) _________________________


FlatField* SteadyStokes :: GetPressure ()
{
    // Returns the pressure field of the receiver.

    return (FlatField*)unknown ;
}


MortaredField* SteadyStokes :: GetVelocity ()
{
    // Returns the velocity field of the receiver.

    return velocity ;
}


#endif

