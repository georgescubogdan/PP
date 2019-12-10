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

#ifndef TIMINTEG
#define TIMINTEG

/*!
 * \file timinteg.hxx
 * \brief Manages time integration formulae
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/garbaged.hxx"
#include "core/param.hxx"
#include "core/vector.hxx"

class SemiDiscreteProblem ;
class Field ;
class FlatField ;

/*! \class TimeIntegrationScheme
   * \brief Manages formulae for time integration of differential problems.
   * 
   * A time-integration scheme is responsible for advancing in time a
   * partial differential problem (the scheme is usually attribute of that
   * problem). When the scheme is created, it is assigned a time-
   * independent problem <{semiDiscreteProblem}>, typically a Helmholtz
   * problem, which it solves at every time step.
   * 
   * The time-related pieces of information are: 
   *   . <{solution}>, the field (flat or mortared) to be advanced in time,
   *   . <{currentTimeStep}>, initially set to 0,
   *   . <{time}>, the current date,
   *   . <{dt}>, the time increment.
   * 
   * In order to save space and time, <{solution}> is shared with the semi-
   * discrete problem (its attribute <{unknown}>).
   * 
   * <{solution}> is a flat or mortared field.
   * By contrast, <{pastSolutions}> contains flat fields, because assembly
   * or mortar operations never need to be applied on them.
   * 
   * The time-integration scheme assigns time-integration related values to
   * <{semiDiscreteProblem}>: a 'lambda' parameter, a flag for implicit or
   * explicit solving mode, a 'weak-history' field and a 'strong-history'
   * field.
   * An additional contribution <{additionalWeakHistory}> may be supplied
   * externally. This is typically the case in Navier-Stokes problems,
   * where the nonlinear term is supplied externally to the scheme
   * responsable of the integration of the diffusive (unsteady Stokes)
   * part.
   * 
   * Some time-integration schemes are not self starting, e.g., BDF2, BDF3,
   * AB2, AB3. The user has two possibilities:
   * - starting with a self-starting scheme (attribute <{starter}>), such
   *   as a BDF1 or Crank-Nicolson (theta method with theta=0.5);
   * - starting with user-supplied values for the starting steps v(-1),
   *   v(-2), etc.
   * Attribute <{useStarter}> is set to "true" (resp. to "false") in the
   * former (resp. latter) case.
   * The type of <{starter}> is defined in the subclasses.
   */ 
class TimeIntegrationScheme : public GarbagedObject
{

protected :
    int                    order;
    int                    currentTimeStep ;
    real                   time ;
    real                   dt ;
    SemiDiscreteProblem*   semiDiscreteProblem ;
    Field*                 solution ;
    Vector<FlatField*>*    pastSolutions ;
    TimeIntegrationScheme* starter ;
    boolean                useStarter ;
    FlatField*             additionalWeakHistory ;

public :

    // = CONSTRUCTORS

    TimeIntegrationScheme () ;
    ~TimeIntegrationScheme () ;

    // = SETTINGS

    void           SetAdditionalWeakHistory (FlatField* hist) ;
    void           SetCurrentTimeStep (int n) ;
    virtual void   SetField (Field* x) ;
    virtual void   SetSemiDiscreteProblem (SemiDiscreteProblem* prob) ;
    void           SetTime (real t) ;
    virtual void   SetTimeIncrement (real deltaT) ;
    real           GetTimeIncrement () ;
    virtual void   SetValues (int i, FlatField* x) ;
    void           UseStarter (boolean b) ;

    // = INQUIRY
  
    real           GetTime () ;
    Field*         GetPastSolution(int q);
    int            GetOrder();

    // = TIME STEPPING

    virtual void   ComputePastSequence (FlatField* y) ;
    virtual void   DoOneStep() = 0 ;
    virtual Field* GetSolution () ;
    virtual void   UpdateTime () = 0 ;


    virtual void PrintPastSequence ();
    // SAVE DATA
//     void      PrintInProc0Past(const char *nom, FlatField* format = NULL);
//     void      ReadFromProc0Past(const char*nom, FlatField* format = NULL) ;
    void      PrintInProc0Past(const char *nom);
    void      ReadFromProc0Past(const char*nom);
} ;


/*! \class AdamsBashforth
   * \brief Manages the Adams-Bashforth time stepping family til order 4.
   * 
   * Integrates approximately
   *   du/dt = F(u,t).
   * At every time step n+1, computes the new value u(n+1). For example, 
   * for a 4-th order scheme:
   *   u(n+1) = u(n) + dt/24 (55 F(n) - 59 F(n-1) + 37 F(n-2) - 9 F(n-3))
   * Attribute <{coefficients}> stores the scalars (55/24, etc).
   * 
   * At any time during any time step n+1 (> 3) the scheme uses 2 +
   * <{order}> field objects:
   *   . <{solution}>, which typically stores u(n+1);
   *   . <{pastSolution}>, which typically stores u(n);
   *   . <{pastFunctionals}>, which stores the <{order}> evaluations of the
   *     functional: F(n), F(n-1), F(n-2), etc.
   * 
   * These same fields are used over and over during the time history, so 
   * as to avoid allocations/deallocations.
   * <{solution}> is not allocated by the scheme: it is received as initial
   * value field; also, it is transmitted at every time step to the semi-
   * discrete problem.
   */ 
class AdamsBashforth : public TimeIntegrationScheme
{

protected :

    RealVector*            coefficients ;
    Vector<FlatField*>*    pastFunctionals ;

    FlatField*             ComputeStrongHistory () ;
    FlatField*             ComputeWeakHistory () ;
    TimeIntegrationScheme* GetStarter () ;
    void                   SetCoefficients () ;
  
public :

    // = CONSTRUCTORS

    AdamsBashforth (int ord) ;
    ~AdamsBashforth () ;
  
    // = SETTINGS
 
    void SetField (Field* x) ;
    void SetSemiDiscreteProblem (SemiDiscreteProblem* prob) ;
    void SetTimeIncrement (real deltaT) ;

    // = TIME STEPPING

    void ComputePastSequence (FlatField* y) ;
    void DoOneStep() ;
    void UpdateTime () ;
} ;


/*! \class BDF
   * \brief Manages the backward-differentiation formula time-stepping family.
   * 
   * Integrates approximately
   *   du/dt = F(u,t).
   * At every time step n+1, computes the new value u(n+1):
   *   a_0/dt u(n+1) = F(n+1) - 1/dt Sum {a_1 u(n) + a_2 u(n-1) + ...}.
   * For example, for BDF2, a_0 = 1.5, a_1 = -2, a_2 = 0.5.
   * 
   * Attribute <{a0}> stores the scalar a_0. <{coefficients}> stores a_1,
   * a_2, a_3.
   * The inherited attribute <{solution}> stores u(n+1).
   * Attribute <{pastSolutions}> contains <{order}> fields: u(n), u(n-1),
   * etc.
   * 
   * A BDF scheme of order > 1 is not self-starting. Then user can provide
   * the BDF with a <{starter}> scheme. By default, see method GetStarter.
   */ 

class BDF : public TimeIntegrationScheme
{

protected :

    real                   a0 ;
    RealVector*            coefficients ;

    FlatField*             ComputeStrongHistory () ;
    TimeIntegrationScheme* GetStarter () ;
    void                   SetCoefficients () ;

public :

    // = CONSTRUCTORS

    BDF (int ord) ;
    ~BDF () ;
  
    // = SETTINGS
 
    void SetTimeIncrement (real deltaT) ;

    // = TIME STEPPING

    void DoOneStep() ;
    void UpdateTime () ;
    void RetrievePastSolutionAtOrder();

    real*  GetCoefficients ();
} ;


/*! \class Extrapolation
   * \brief Manages extrapolation schemes.
   * 
   * Currently, 1-point to 4-point extrapolation is implemented.
   * For example, the 3-point scheme implements the formula
   *     x(n+1) = 3 x(n) - 3 x(n-1) + x(n-2).
   * 
   * This is not truly a time-integration scheme. It is just an 
   * extrapolation formula, implemented as a time-integration scheme for
   * convenience. Most of the attributes and methods inherited from
   * TimeIntegrationScheme are not used.
   */ 

class Extrapolation : public TimeIntegrationScheme
{

protected :

    int         nbPoints ;
    RealVector* coefficients ;

    void SetCoefficients () ;

public :

    // = CONSTRUCTORS

    Extrapolation (int np) ;
    ~Extrapolation () ;
  
    // = SETTINGS
 
    void SetField (Field* x) ;
    void SetSemiDiscreteProblem (SemiDiscreteProblem* prob) ;

    // = TIME STEPPING

    void ComputePastSequence (FlatField* y) ;
    void SetValues (int i, FlatField* x) ;
    void DoOneStep() ;
    void UpdateTime () ;
    void PrintPastSequence() ;
} ;



/*! \class RungeKutta
   * \brief Manages the Runge-Kutta time stepping schemes of order 2 and 4.
   * 
   * Integrates approximately:
   *   du/dt = F(u,t).
   * At every time step n+1, computes the new value u(n+1). For example, 
   * for a 4-th order scheme:
   *   u(n+1) = u(n) + dt/6 (k1 + 2 k2 + 2 k3 + k4)
   * 
   * At any time during any time step n+1 the scheme uses 2 +
   * <{order}> field objects:
   *   . <{solution}>, which typically stores u(n+1);
   *   . <{k1}>, ..., for the k_i fields;
   *   . <{pastSolution}>, which typically stores u(n);
   * These same fields are used over and over during the time history, so 
   * as to avoid allocations/deallocations.
   * <{solution}> is not allocated by the scheme: it is received as initial
   * value field; also, it is transmitted at every time step to the semi-
   * discrete problem.
   * 
   * Space occupation can be improved: <{uwork}> can be avoided and the
   * number of k_i's reduced.
   * 
   * Warning. Do not try and use a Runge-Kutta scheme for a (Navier-)
   * Stokes problem. The calculation of k1, k2, etc, depends on intermiate-
   * step values of the pressure field which are not implemented.
   */ 

class RungeKutta : public TimeIntegrationScheme
{

protected :

    FlatField* k1 ;
    FlatField* k2 ;
    FlatField* k3 ;
    FlatField* k4 ;
    FlatField* uwork ;

    FlatField* ComputeStrongHistory () ;
    FlatField* ComputeWeakHistory () ;

public :

    // = CONSTRUCTORS

    RungeKutta (int ord) ;
    ~RungeKutta () ;
  
    // = SETTINGS
 
    void SetField (Field* x) ;
    void SetSemiDiscreteProblem (SemiDiscreteProblem* prob) ;

    // = TIME STEPPING

    void DoOneStep() ;
    void UpdateTime () ;
} ;


/*! \class ThetaMethod
   * \brief Manages the generalized trapezoidal rule family.
   * 
   * Integrates approximately:
   *   du/dt = F(u,t).
   * by means of
   *   u(n+1) = u(n) + dt { (1-theta) F(u(n)) + theta F(u(n+1)) }.
   * For example, values for <{theta}> equal to 0, 1/2 and 1 yield the 
   * Euler forward (explicit), Crank-Nicholson (trapezoidal rule) and 
   * backward Euler algorithms, respectively.
   * 
   * In order to save memory space, <{F}> stores alternatively F(u(n)) and
   * F(u(n+1)).
   */

class ThetaMethod : public TimeIntegrationScheme
{

protected :

    real       theta ;
    FlatField* F ;

    FlatField* ComputeStrongHistory () ;
    FlatField* ComputeWeakHistory () ;

public :

    // = CONSTRUCTORS

    ThetaMethod (real thet) ;
    ~ThetaMethod () ;
  
    // = SETTINGS
 
    void SetField (Field* x) ;
    void SetSemiDiscreteProblem (SemiDiscreteProblem* prob) ;

    // = TIME STEPPING

    void DoOneStep() ;
    void UpdateTime () ;
} ;


#endif

