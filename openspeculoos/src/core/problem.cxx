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

// problem.cxx

#include "core/problem.hxx"
#include "core/mesh.hxx"
#include "core/field.hxx"
#include "core/timinteg.hxx"
#include "core/vector.hxx"
#include "core/point.hxx"
#include "core/parent.hxx"
#include "core/matrix.hxx"
#include "core/bc.hxx"
#include "core/solver.hxx"
#include "core/precond.hxx"
#include "core/funct.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include <math.h>
#include <stdio.h>



// Ajout Vincent Keller (needs libmonitor)
//#include "../MPIMonitor_v1.1/mpimonitor.h"



//________________________________ Problem ____________________________________


Problem :: Problem (Mesh* m) 
    : GarbagedObject ()
{
    // Constructor. Initializes the receiver as a problem defined on 'm'.
  
    ref(m) ;
    mesh                = m ;
    dirichletConditions = new Vector<DirichletCondition*>() ;
    neumannConditions   = new Vector<NeumannCondition*>() ;
    forcing             = NULL ;
    forcingFunction     = NULL ;
    time                = ZERO ;
}


Problem :: ~Problem ()
{
    // Destructor.

    DeleteDirichletConditions() ;
    DeleteNeumannConditions() ;

    unref(mesh) ;
    unref(forcing) ;
    unref(forcingFunction) ;
    //  fclose(ficdeb);
}


void Problem :: AddDirichletCondition (DirichletCondition* dc)
{
    // Adds 'dc' to the list of Dirichlet conditions of the receiver.

    ref(dc) ;
    dirichletConditions->Put(dc) ;
}


void Problem :: AddNeumannCondition (NeumannCondition* nc)
{
    // Adds 'dc' to the list of Neumann conditions of the receiver.

    ref(nc) ;
    neumannConditions->Put(nc) ;
}


void Problem :: ApplyDirichletConditions (Field* x)
{
    // Applies the Dirichlet conditions of the receiver to the relevant component
    // of 'x'.

    int i, n ;

    n = dirichletConditions->GetSize() ;
    for (i=1 ; i <= n ; i++)
	dirichletConditions->At(i)->ApplyOn(x,time) ;
}


void Problem :: ApplyHomogeneousDirichletConditions (Field* x)
{
    // Applies the Dirichlet conditions of the receiver to 'x'. This is done in a
    // homogeneous way (the prescribed values are replaced by 0).
    // The choice of the affected component of 'x' is governed by the Dirichlet
    // conditions.

    int i, n ;

    n = dirichletConditions->GetSize() ;
    for (i=1 ; i <= n ; i++)
	dirichletConditions->At(i)->ApplyHomogeneousOn(x) ;
}


void Problem :: ApplyHomogeneousDirichletConditions (Field* x, int icx)
{
    // Applies the Dirichlet conditions of the receiver to the 'icx'-th component
    // of 'x'. This is done in a homogeneous way (the prescribed values are 
    // replaced by 0).

    int i, n ;

    n = dirichletConditions->GetSize() ;
    for (i=1 ; i <= n ; i++)
	dirichletConditions->At(i)->ApplyHomogeneousOn(x,icx) ;
}


void Problem :: DeleteDirichletConditions ()
{
    // Removes the Dirichlet conditions of the receiver.

    int i ;

    if (dirichletConditions) {
	for (i=1 ; i <= dirichletConditions->GetSize() ; i++)
	    unref(dirichletConditions->At(i)) ;
	delete dirichletConditions ;}
    dirichletConditions = new Vector<DirichletCondition*>() ;
}


void Problem :: DeleteNeumannConditions ()
{
    // Removes the Neumann conditions of the receiver.

    int i ;

    if (neumannConditions) {
	for (i=1 ; i <= neumannConditions->GetSize() ; i++)
	    unref(neumannConditions->At(i)) ;
	delete neumannConditions ;}
    neumannConditions = new Vector<NeumannCondition*>() ;
}


void Problem :: EvaluateForcing (real t)
{
    // Initializes the forcing-term field f to its value at time 't'.
    // Does nothing if the receiver has no function for the forcing-term field.

# ifdef REQUIRE
    Require("forcing term exists", forcing != NULL) ;
    InMethod("Problem::EvaluateForcing(t)") ;
# endif

    if (forcingFunction)
	forcing->SetTo(forcingFunction,t) ;
}


FlatField* Problem :: GetForcing ()
{
    // Returns the forcing-term of the receiver.

    return forcing ;
}


Vector<DirichletCondition*>* Problem :: GetDirichletConditions ()
{
    // Returns the vector of the Dirichlet conditions of the receiver.

    return dirichletConditions ;
}


Vector<DirichletCondition*>* Problem :: GetHomogeneousDirichletConditions ()
{
    // Returns a new vector containing copies of the Dirichlet conditions with
    // all coefficients of the prescribed-value fields set to 0.

    Vector<DirichletCondition*> *answer ;
    int                         i ;

    answer = new Vector<DirichletCondition*>() ;

    for (i=1 ; i <= dirichletConditions->GetSize() ; i++)
	answer->Put(dirichletConditions->At(i)->DuplicateHomogeneous()) ;

    return answer ;
}
    
Vector<NeumannCondition*>*  Problem :: GetNeumannConditions()
{
    // Returns the vector of the Neumann conditions of the receiver.
    return neumannConditions;
}
    

Mesh* Problem :: GetMesh ()
{
    // Returns the mesh of the receiver.

    return mesh ;
}


void Problem :: SetDirichletConditions (Vector<DirichletCondition*>* dcs)
{
    // Sets the Dirichlet conditions of the receiver to 'dcs'.

    DirichletCondition *dc ;
    int                i ;

    DeleteDirichletConditions() ;

    for (i=1 ; i <= dcs->GetSize() ; i++) {
	dc = dcs->At(i) ;
	ref(dc) ;
	dirichletConditions->Put(dc) ;    
    }
}


void Problem :: SetForcing (FlatField* f)
{
    // Sets the forcing term of the receiver to 'f'.
    // 'f' may be NULL.

# ifdef REQUIRE
    Require("same mesh", f == NULL || f->GetMesh() == mesh) ;
    InMethod("Problem::SetForcing(f)") ;
# endif

    ref(f) ;
    unref(forcing) ;
    forcing = f ;
}


void Problem :: SetForcingFunction (Function* funct)
{
    // Sets the function for the forcing term of the receiver to 'funct'.
    // 'funct' may be NULL.

    ref(funct) ;
    unref(forcingFunction) ;
    forcingFunction = funct ;
}


void Problem :: SetNeumannConditions (Vector<NeumannCondition*>* ncs)
{
    // Sets the Neumann conditions of the receiver to 'ncs'.

    NeumannCondition *nc ;
    int              i ;

    DeleteNeumannConditions() ;
    for (i=1 ; i <= ncs->GetSize() ; i++) {
	nc = ncs->At(i) ;
	ref(nc) ;
	neumannConditions->Put(nc) ;    
    }
}


void Problem :: SetTime (real t)
{
    // Sets the current time to 't'.

    time = t ;
}


//_______________________________ HeatEquation ________________________________


HeatEquation :: HeatEquation (Mesh* m)
    : Problem(m)
{
    // Constructor. Creates a heat-equation problem on 'm'.

    u                     = NULL ;
    timeIntegrationScheme = NULL ;
    helmholtz             = new Helmholtz(m) ;
}


HeatEquation :: ~HeatEquation ()
{
    // Destructor.

    unref(u) ;
    unref(timeIntegrationScheme) ;
    unref(helmholtz) ;
}


void HeatEquation :: AddDirichletCondition (DirichletCondition* dc)
{
    // Adds 'dc' to the list of Dirichlet conditions of the receiver.

    helmholtz->AddDirichletCondition(dc) ;
}


void HeatEquation :: AddNeumannCondition (NeumannCondition* nc)
{
    // Adds 'dc' to the list of Neumann conditions of the receiver.

    helmholtz->AddNeumannCondition(nc) ;
}


real Helmholtz :: GetC1 ()
{
    // Returns the c1 coefficient of the Couzy preconditioner
    //   inv(P) = c1 inv(B) + c2 inv(A).

    return lambda ;
}


real Helmholtz :: GetC2 ()
{
    // Returns the c2 coefficient of the Couzy preconditioner
    //   inv(P) = c1 inv(B) + c2 inv(A).

    return alpha ;
}


Helmholtz* HeatEquation :: GetHelmholtz ()
{
    // A debugging tool.
    // Returns the Helmholtz subproblem of the receiver.

    return helmholtz ;
}


Field* HeatEquation :: GetUnknown ()
{
    // Returns the current value of the unknown field.

    return u ;
}


void HeatEquation :: SetForcing (FlatField* f)
{
    // Sets the forcing term of the receiver to 'f'.
    // 'f' may be NULL.

# ifdef REQUIRE
    Require("same mesh", f == NULL || f->GetMesh() == mesh) ;
    InMethod("HeatEquation::SetForcing(f)") ;
# endif

    helmholtz->SetForcing(f) ;
}


void HeatEquation :: SetForcingFunction (Function* funct)
{
    // Sets the function for the forcing term of the receiver to 'funct'.
    // 'funct' may be NULL.

    helmholtz->SetForcingFunction(funct) ;
}


void HeatEquation :: SetField (Field* x)
{
    // Assigns 'x' as unknown field for the receiver.

# ifdef REQUIRE
    Require("same mesh", x->GetMesh() == mesh) ;
    InMethod("HeatEquation::SetField(x)") ;
# endif

    ref(x) ;
    unref(u) ;
    u = x ;

    if (timeIntegrationScheme)
	timeIntegrationScheme->SetField(x) ;
}


void HeatEquation :: SetIntegrationRule (FlatField* rule)
{
    // Sets the integration rule of the receiver to 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("HeatEquation::SetIntegrationRule(rule)") ;
# endif

    helmholtz->SetIntegrationRule(rule) ;
}


void HeatEquation :: SetPreconditioner (Preconditioner* p)
{
    // Assigns 'p' as preconditioner for any equation system to be solved. 'p'
    // may be NULL.

    helmholtz->SetPreconditioner(p) ;
}


void HeatEquation :: SetSolver (Solver* s)
{
    // Assigns 's' as solver for any equation system to be solved.

    helmholtz->SetSolver(s) ;
}


void HeatEquation :: SetThermalDiffusivity (real val)
{
    // Sets the C*C coefficient of the receiver to 'val'.

    helmholtz->SetAlpha(val) ;
}


void HeatEquation :: SetTime (real t)
{
    // Sets the clock of the receiver to 't'. Used when evaluating the forcing
    // term and the time-dependent Dirichlet conditions.

    time = t ;

    helmholtz->SetTime(t) ;
}


void HeatEquation :: SetTimeIntegrationScheme (TimeIntegrationScheme* scheme)
{
    // Assigns 'scheme' as time-integration scheme for the receiver. The receiver
    // in turn assigns its Helmholtz problem to 'scheme'.

    ref(scheme) ;
    unref(timeIntegrationScheme) ;
    timeIntegrationScheme = scheme ;

    scheme->SetSemiDiscreteProblem(helmholtz) ;
    if (u)
	scheme->SetField(u) ;
}


void HeatEquation :: Solve (int nstep)
{
    // Solves 'nstep' time steps.

# ifdef REQUIRE
    Require("has a time-integration scheme", timeIntegrationScheme != NULL) ;
    Require("has initial values", u != NULL) ;
    InMethod("HeatEquation::Solve(nstep)") ;
# endif

    int i ;

    for (i=1 ; i<=nstep ; i++)
	timeIntegrationScheme->DoOneStep() ;
}


//__________________________________ Stokes ___________________________________


Stokes :: Stokes (Mesh* m, DecompositionType typ)
    : Problem(m)
{
    // Constructor. Creates an unsteady Stokes problem on 'm'. 
    // The type of the generalized-LU decomposition of the underlying steady 
    // Stokes problem is given by 'typ'.

# ifdef REQUIRE
    Require("valid decomposition type", typ==UZAWA || typ==BP1 || typ==BP3 ||
	    typ==BP1_PC || typ==BP3_PC || typ==YOSIDA) ;
    InMethod("Stokes::Stokes(m,typ)") ;
# endif

    velocity              = NULL ;
    pressure              = NULL ;
    timeIntegrationScheme = NULL ;
    steadyStokes          = new SteadyStokes(m,typ) ;
}


Stokes :: ~Stokes ()
{
    // Destructor.

    unref(timeIntegrationScheme) ;
    unref(steadyStokes) ;
    unref(velocity) ;
    unref(pressure) ;
}


void Stokes :: AddDirichletConditionV (DirichletCondition* dc)
{
    // Adds 'dc' to the list of Dirichlet conditions on the velocity.

    steadyStokes->AddDirichletConditionV(dc) ;
}


void Stokes :: AddNeumannConditionV (NeumannCondition* nc)
{
    // Adds 'nc' to the list of Neumann conditions on the velocity.

    steadyStokes->AddNeumannConditionV(nc) ;
}


SteadyStokes* Stokes :: GetSteadyStokes ()
{
    // Returns the steady Stokes subproblem of the receiver.

    return steadyStokes ;
}


FlatField* Stokes :: GetPressure ()
{
    // Returns the current value of the pressure field.

    return pressure ;
}


MortaredField* Stokes :: GetVelocity ()
{
    // Returns the current value of the velocity field.

    return velocity ;
}


void Stokes :: SetForcing (FlatField* f)
{
    // Sets the forcing term of the receiver to 'f'.
    // 'f' may be NULL.

# ifdef REQUIRE
    Require("same mesh", f == NULL || f->GetMesh() == mesh) ;
    InMethod("Stokes::SetForcing(f)") ;
# endif

    steadyStokes->SetForcing(f) ;
}


void Stokes :: SetForcingFunction (Function* funct)
{
    // Sets the function for the forcing term of the receiver to 'funct'.
    // 'funct' may be NULL.

    steadyStokes->SetForcingFunction(funct) ;
}


void Stokes :: SetIntegrationRuleMass (FlatField* rule)
{
    // Initializes the integration rule for the mass equation to 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("Stokes::SetIntegrationRuleMass(rule)") ;
# endif

    steadyStokes->SetIntegrationRuleMass(rule) ;
}


void Stokes :: SetIntegrationRuleMomentum (FlatField* rule)
{
    // Initializes the integration rule for the momentum equations to 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("Stokes::SetIntegrationRuleMomentum(rule)") ;
# endif

    steadyStokes->SetIntegrationRuleMomentum(rule) ;
}


void Stokes :: SetPressure (FlatField* p)
{
    // Assigns 'p' as pressure field for the receiver.
    // Only the interpolation of 'p' is relevant; the values may contain garbage.
    // This very object 'p' will contain the current values of the pressure field
    // at any time throughout the analysis, so do not delete it!

# ifdef REQUIRE
    Require("same mesh", p->GetMesh() == mesh) ;
    InMethod("Stokes::SetPressure(p)") ;
# endif

    ref(p) ;
    unref(pressure) ;
    pressure = p ;

    steadyStokes->SetPressure(p) ;
}


void Stokes :: SetVelocity (MortaredField* v)
{
    // Assigns 'v' as velocity field for the receiver.
    // Only the interpolation of 'v' is relevant; the values may contain garbage.
    // This very object 'v' will contain the current values of the velocity field
    // at any time throughout the analysis, so do not delete it!

# ifdef REQUIRE
    Require("same mesh", v->GetMesh() == mesh) ;
    InMethod("Stokes::SetVelocity(v)") ;
# endif

    ref(v) ;
    unref(velocity) ;
    velocity = v ;

    steadyStokes->SetVelocity(v) ;

    if (timeIntegrationScheme)
	timeIntegrationScheme->SetField(v) ;
}


void Stokes :: SetSolver (Solver* s)
{
    // Assigns 's' as solver for any equation system to be solved.

    steadyStokes->SetSolver(s) ;
}


void Stokes :: SetNu (real val)
{
    // Sets the kinematic viscosity of the receiver to 'val'.

    steadyStokes->SetNu(val) ;
}

void Stokes :: ImposeNoPressureNullSpace()
{
    steadyStokes->ImposeNoPressureNullSpace();
}

void Stokes :: SetTimeIntegrationScheme (TimeIntegrationScheme* scheme)
{
    // Assigns 'scheme' as time-integration scheme for the receiver. The receiver
    // in turn assigns its steady Stokes problem to 'scheme'.

    ref(scheme) ;
    unref(timeIntegrationScheme) ;
    timeIntegrationScheme = scheme ;

    scheme->SetSemiDiscreteProblem(steadyStokes) ;
    if (velocity)
	scheme->SetField(velocity) ;
}

void Stokes :: SetValues (int i, FlatField* x)
{
    // Sets to 'x' the starting values of the velocity at step '-i'.
    // For example, for a BDF3/EX3 integration, i=0,1,2 initializes respectively
    // u(0), u(-1) and u(-2).

# ifdef REQUIRE
    Require("interpolation of velocity already set", velocity != NULL) ;
    Require("has diffusion scheme", timeIntegrationScheme != NULL) ;
    Require("has integration rule", steadyStokes->GetIntegrationRuleMomentum()
	    != NULL) ;
    InMethod("NavierStokes::SetPastValues(i,x)") ;
# endif

    // Stokes part (stores the velocity)
    timeIntegrationScheme->SetValues(i,x) ;
}


void Stokes :: Solve (int nstep)
{
    // Solves 'nstep' time steps.

# ifdef REQUIRE
    Require("has a time-integration scheme", timeIntegrationScheme != NULL) ;
    Require("has initial velocity", velocity != NULL) ;
    Require("has interpolation for initial pressure", pressure != NULL) ;
    InMethod("Stokes::Solve(nstep)") ;
# endif

    int i ;

    for (i=1 ; i<=nstep ; i++) {
	//    printf("\n\n Start time step %d\n",i) ;
	timeIntegrationScheme->DoOneStep() ;
    }
}

//________________________________ NavierStokes _______________________________


NavierStokes :: NavierStokes (Mesh* m, DecompositionType typ)
    : Stokes(m,typ)
{
    // Constructor. Creates a Navier-Stokes problem on 'm'.
    // The type of the generalized-LU decomposition of the underlying steady 
    // Stokes problem is given by 'typ'.

    convectionScheme = NULL ;
    nonlin           = NULL ;
}


NavierStokes :: ~NavierStokes ()
{
    // Destructor.

    unref(convectionScheme) ;
    unref(nonlin) ;
}


void NavierStokes :: GetNonlinearTerm (FlatField* y)
{
    // Updates the time to step n+1, then initializes 'y' to the weak form of the
    // extrapolated nonlinear term NL(n+1).
    //
    // Involves the calculation of the exact nonlinear term at the previous step:
    //   NL(n) = -v(n).grad v(n).
  
# ifdef REQUIRE
    Require("has integration rule", steadyStokes->GetIntegrationRuleMomentum()
	    != NULL) ;
    InMethod("NavierStokes::GetNonlinearTerm(y)") ;
# endif

    FlatField *v ;

    // compute -v.grad v at last time step n (uses y as temporary storage)
    v = velocity->GetMain() ;
    y->SetToWeakConvection(v,v,steadyStokes->GetIntegrationRuleMomentum()) ;
    y->SwitchSigns() ;
    convectionScheme->UpdateTime() ;
    convectionScheme->SetValues(1,y) ;

    // compute -v.grad v at new time step n+1
    convectionScheme->ComputePastSequence(y) ;   // y = approximation of NL(n+1)
}


void NavierStokes :: SetConvectionScheme (TimeIntegrationScheme* scheme)
{
    // Assigns 'scheme' as time-integration scheme for the convective part of the
    // receiver.
    // Presently 'scheme' must be an Extrapolation (note: this is not checked).

    ref(scheme) ;
    unref(convectionScheme) ;
    convectionScheme = scheme ;

    if (velocity)
	scheme->SetField(velocity) ;
}


void NavierStokes :: SetReynolds (real re)
{
    // Initializes the Reynolds number of the receiver to 're'.

# ifdef REQUIRE
    Require("valid Reynolds number", re > 0.) ;
    InMethod("NavierStokes::SetReynolds(re)") ;
# endif

    SetNu(ONE/re) ;
}


void NavierStokes :: SetValues (int i, FlatField* x)
{
    // Sets to 'x' the starting values of the velocity at step '-i'.
    // For example, for a BDF3/EX3 integration, i=0,1,2 initializes respectively
    // u(0), u(-1) and u(-2).

# ifdef REQUIRE
    Require("interpolation of velocity already set", velocity != NULL) ;
    Require("has diffusion scheme", timeIntegrationScheme != NULL) ;
    Require("has convection scheme", convectionScheme != NULL) ;
    Require("has integration rule", steadyStokes->GetIntegrationRuleMomentum()
	    != NULL) ;
    InMethod("NavierStokes::SetPastValues(i,x)") ;
# endif

    // Stokes part (stores the velocity)
    timeIntegrationScheme->SetValues(i,x) ;

    // convection part (stores -v.grad v)
    if (i != 0) {
	nonlin->SetToWeakConvection(x,x, steadyStokes->GetIntegrationRuleMomentum()) ;
	nonlin->SwitchSigns() ;
	convectionScheme->SetValues(i,nonlin) ;
    }
}


void NavierStokes :: SetVelocity (MortaredField* v)
{
    // Assigns 'v' as velocity field for the receiver.
    // Only the interpolation of 'v' is relevant; the values may contain garbage.
    // This very object 'v' will contain the current values of the velocity field
    // at any time throughout the analysis, so do not delete it!

# ifdef REQUIRE
    Require("same mesh", v->GetMesh() == mesh) ;
    InMethod("Stokes::SetVelocity(v)") ;
# endif

    Stokes::SetVelocity(v) ;

    if (convectionScheme)
	convectionScheme->SetField(v) ;

    unref(nonlin) ;
    nonlin = v->GetMain()->DuplicateEmpty() ;
}


//void NavierStokes :: Solve (int nstep)
real NavierStokes :: Solve (int nstep)
{
    // Solves 'nstep' time steps.

# ifdef REQUIRE
    Require("has diffusion scheme", timeIntegrationScheme != NULL) ;
    Require("has convection scheme", convectionScheme != NULL) ;
    Require("has interpolation for velocity", velocity != NULL) ;
    Require("has interpolation for pressure", pressure != NULL) ;
    InMethod("NavierStokes::Solve(nstep)") ;
# endif

    int i ;

    for (i=1 ; i<=nstep ; i++) {
	// convective term
	GetNonlinearTerm(nonlin) ;
	timeIntegrationScheme->SetAdditionalWeakHistory(nonlin) ;
    
	// Stokes part
	timeIntegrationScheme->DoOneStep() ;
    }

    // compute the infinite norm of v_n+1 - v_n
    real vnp1mvn = 0.;
    FlatField* tempo, *vhat;

    vhat = (FlatField *)timeIntegrationScheme->GetPastSolution(1);
    tempo = vhat->Duplicate();
    tempo->SwitchSigns();
    tempo->Add(velocity->GetMain());
    vnp1mvn = tempo->GetInfiniteNorm();

    delete tempo;

    return vnp1mvn;
}


//_____________________________ SemiDiscreteProblem ___________________________


SemiDiscreteProblem :: SemiDiscreteProblem (Mesh* m)
    : Problem (m)
{
    // Constructor. Initializes the receiver as a problem defined on 'm'.

    lambda          = ZERO ;
    unknown         = NULL ;
    integrationRule = NULL ;
    solver          = NULL ;
    weakHistory     = NULL ;
    strongHistory   = NULL ;
    mode            = IMPLICIT ;
    preconditioner  = NULL ;
    inverseB        = NULL ;
}


SemiDiscreteProblem :: ~SemiDiscreteProblem ()
{
    // Destructor.
 
    unref(weakHistory) ;
    unref(strongHistory) ;
    unref(solver) ;
    unref(preconditioner) ;
    unref(inverseB) ;
    unref(unknown) ;
    unref(integrationRule) ;
}


void SemiDiscreteProblem :: ApplyInverseB (Field* y) 
{
    // Initializes 'y' to 'inv(B) y', where B is the mass (i.e., discrete weak 
    // identity) operator.
    // If 'y' and the unknown field are mortared fields, inv(B) encompasses 
    // mortar constraints and homogeneous Dirichlet conditions.

# ifdef REQUIRE
    Require("'y' has same type as inv(B)", y->HasSameTypeAs(unknown)) ;
    Require("'y' has same interpolation as inv(B)", 
	    y->HasSameInterpolationAs(unknown)) ;
    InMethod("SemiDiscreteProblem::ApplyInverseB(y) ") ;
# endif

    int i ;

    // 1. Assemble y

    y->Assemble() ;
  
    // 2. Initialize y to inv(B) y

    GetInverseB() ;
    for (i=1 ; i <= y->GetNbComponents() ; i++)
	y->Multiply(inverseB,1,i) ;

    // 3. De-assemble y

    y->Distribute() ;
} 


void SemiDiscreteProblem :: DeleteInverseB ()
{
    // Deletes inv(B), if already computed.
    // Useful for "synchronizing" the interpolation of inv(B) with that of the
    // unknown field, if the latter has been modified. 
    // Can be safely applied at any time: if inv(B) is requested later on, it 
    // will be automatically recalculated.

    unref(inverseB) ;
    inverseB = NULL ;
}


real SemiDiscreteProblem :: GetC1 ()
{
    // Issues an error message.

    NotImplemented("SemiDiscreteProblem::GetC1()") ;
    return 0. ;
}


real SemiDiscreteProblem :: GetC2 ()
{
    // Issues an error message.

    NotImplemented("SemiDiscreteProblem::GetC2()") ;
    return 0. ;
}


FlatField* SemiDiscreteProblem :: GetIntegrationRule ()
{
    // Returns the integration rule of the receiver.

    return integrationRule ;
}


Field* SemiDiscreteProblem :: GetInverseB ()
{
    // Returns inv(B), the inverse of the assembled mass matrix. Since the mass 
    // matrix is diagonal, it is stored as a field.
    //
    // inv(B) has the same nature (flat or mortared) and same interpolation as
    // the unknown field. It has 1 component.
    //
    // If inv(B) is a mortared field, only its values at the internal dofs are
    // meaningful, to the assembly operation. So do not multiply it with a flat 
    // field!
    //
    // The homogeneous Dirichlet conditions are accounted for.

# ifdef REQUIRE
    Require("the unknown field exists", unknown != NULL) ;
    Require("same interpolation", inverseB == NULL || 
	    inverseB->HasSameInterpolationAs(unknown)) ;
    InMethod("SemiDiscreteProblem::GetInverseB()") ;
# endif

    if (inverseB == NULL) {
	inverseB = unknown->FieldDuplicateEmpty(1) ;              // compute B
	inverseB->GetMain()->CopyInterpolateFrom(mesh->GetJacobian()) ;
	inverseB->GetMain()->MultiplyByWeights() ;
	inverseB->Assemble() ;
	inverseB->Inverse() ;                                     // compute inv(B)
	ApplyHomogeneousDirichletConditions(inverseB,1) ;
    }

    return inverseB ;
}


real SemiDiscreteProblem :: GetLambda ()
{
    // Returns the lambda coefficient of the receiver.
 
    return lambda ;
}


FullMatrix* SemiDiscreteProblem :: GetMatrixOperator ()
{
    // Returns the matrix operator A of the receiver. Constructs it by 'probing'
    // (expensive!).
    // The order of A is equal to the number of degrees of freedom of the mesh.
    // Boundary conditions and mortar constraints are taken into account.

    Field      *temp1, *temp2 ;
    FullMatrix *A ;
    int        i, nbdof ;

    // 1. Compute A


    nbdof = unknown->GetNbDofs() ;

    temp1 = unknown->FieldDuplicateEmpty() ;
    temp2 = unknown->FieldDuplicateEmpty() ;

    printf(" Nbdofs = %d \n",nbdof);

    A     = new FullMatrix(nbdof,nbdof) ;

    for (i=1 ; i<=nbdof ; i++) {
	if (i>1) 
	    temp1->SetDof(i-1,ZERO) ;
	temp1->SetDof(i,ONE) ;
	ApplyHomogeneousDirichletConditions(temp1) ;
	temp1->Distribute() ;
	GetAx(temp1,temp2);
	temp2->Assemble() ;
	ApplyHomogeneousDirichletConditions(temp2) ;
	temp2->CopyDofsToMatrixColumn(A,i) ;
    }

    // 2. Check symmetry of A

    // A->PrintMatlab() ; 

    if (A->TestSymmetry(true))
	printf("\n  The matrix operator is symmetric\n") ;
    else {
	printf("\n  The matrix operator is NOT symmetric\n") ;
	Exit() ;
    }

    // 3. Take homogeneous Dirichlet bc's into account in LHS, using a trick for
    //    detecting where Dirichlet bc's apply

    temp1->SetValues(ONE) ;
    ApplyHomogeneousDirichletConditions(temp1) ;

    for (i=1 ; i<=nbdof ; i++)
	if (temp1->GetDof(i) != ONE)
	    A->At(i,i) = ONE ;

    delete temp1 ;
    delete temp2 ;

    return A ;
}


Preconditioner* SemiDiscreteProblem :: GetPreconditioner ()
{
    // Returns the preconditioner of the receiver.

    return preconditioner ;
}


Solver* SemiDiscreteProblem :: GetSolver ()
{
    // Returns the equation solver of the receiver. Creates it if its does not
    // exist yet.

    if (solver == NULL) {
	solver = new ConjugateGradient(500,DEFAULT_TOLERANCE) ;
	solver->SetVerbosity(3) ;
    }
  
    return solver ;
}


Field* SemiDiscreteProblem :: GetUnknown ()
{
    // Returns the field of unknowns of the receiver.

    return unknown ;
}


boolean SemiDiscreteProblem :: IsHelmholtz () 
{
    // Return 'false', since by default the receiver's type is not Helmholtz.

    return false ;
}


boolean SemiDiscreteProblem :: IsSteadyStokes () 
{
    // Return 'false', since by default the receiver's type is not SteadyStokes.

    return false ;
}


void SemiDiscreteProblem :: SetIntegrationRule (FlatField* rule)
{
    // Sets the integration rule of the receiver 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("SemiDiscreteProblem::SetIntegrationRule(rule)") ;
# endif

    ref(rule) ;
    unref(integrationRule) ;

    integrationRule = rule ;
}


void SemiDiscreteProblem :: SetLambda (real val)
{
    // Sets the lambda coefficient of the receiver to 'val'.

    lambda = val ;
}


void SemiDiscreteProblem :: SetPreconditioner (Preconditioner* p)
{
    // Sets to 'p' the preconditioner of the receiver. 'p' may be NULL.

    unref(preconditioner) ;
    ref(p) ;

    if (p)
	p->SetProblem(this) ;

    preconditioner = p ;
}


void SemiDiscreteProblem :: SetSolver (Solver* s)
{
    // Sets the equation solver of the receiver to 's'.
    // 's' may be NULL, in which case the receiver will use its default solver.

    unref(solver) ;
    ref(s) ;
    solver = s ;
}


void SemiDiscreteProblem :: SetStrongHistory (FlatField* hist)
{
    // Sets the strong part (the one that does need to be integrated) of the 
    // history field of the receiver to 'hist'.
    // 'hist' can be NULL.
    //
    // Remark. In order to speed up calculations and to save memory space, no 
    //         copy of 'hist' is made, so do not modify 'hist' before the 
    //         receiver is solved.

# ifdef REQUIRE
    Require("same mesh", hist == NULL || hist->GetMesh() == mesh) ;
    Require("same interpolation", hist == NULL || 
	    hist->HasSameInterpolationAs(unknown->GetMain())) ;
    InMethod("SemiDiscreteProblem::SetStrongHistoryField(hist)") ;
# endif

    ref(hist) ;
    unref(strongHistory) ;
    strongHistory = hist ;
}


void SemiDiscreteProblem :: SetTemporalUnknown (Field* u)
{
    // Default implementation: sets the unknown field of the receiver to 'u'

# ifdef REQUIRE
    Require("same mesh", u->GetMesh() == mesh) ;
    InMethod("SemiDiscreteProblem::SetTemporalUnknown(u)") ;
# endif

    SetUnknown(u) ;
}


void SemiDiscreteProblem :: SetToExplicit ()
{
    // Whenever the receiver's operator is to calculated, ignores the functional
    // F(u).
    // Useful if the receiver is the semi-discrete problem associated with a PDE
    // problem solved by an explicit time-integration scheme.

    mode = EXPLICIT ;
}


void SemiDiscreteProblem :: SetToImplicit ()
{
    // Whenever the receiver's operator is to calculated, considers the 
    // functional F(u). This is the default mode.

    mode = IMPLICIT ;
}


void SemiDiscreteProblem :: SetUnknown (Field* u)
{
    // Initializes the unknown field of the receiver to 'u'.
    // 'u' may be NULL.

# ifdef REQUIRE
    Require("same mesh", u == NULL || u->GetMesh() == mesh) ;
    Require("same interpolation", unknown == NULL || u == NULL ||
	    u->HasSameInterpolationAs(unknown)) ;
    InMethod("SemiDiscreteProblem::SetUnknown(u)") ;
# endif

    ref(u) ;
    unref(unknown) ;
    unknown = u ;
}


void SemiDiscreteProblem :: SetWeakHistory (FlatField* hist)
{
    // Sets the weak part (the one that does not need to be integrated) of the 
    // history field of the receiver to 'hist'.
    // 'hist' can be NULL.
    //
    // Remark. In order to speed up calculations and to save memory space, no 
    //         copy of 'hist' is made, so do not modify 'hist' before the 
    //         receiver is solved.

# ifdef REQUIRE
    Require("same mesh", hist == NULL || hist->GetMesh() == mesh) ;
    Require("same interpolation", hist == NULL || 
	    hist->HasSameInterpolationAs(unknown->GetMain())) ;
    InMethod("SemiDiscreteProblem::SetWeakHistory(hist)") ;
# endif

    ref(hist) ;
    unref(weakHistory) ;
    weakHistory = hist ;
}


void SemiDiscreteProblem :: Solve ()
{
    // Solves the receiver.
    // 1. solve

    GetSolver()->Solve(this) ;

    // 2. account for Dirichlet conditions

    ApplyDirichletConditions(unknown) ;

    // 3. de-assemble the unknown field

    unknown->Distribute() ;
}


//_________________________________ Helmholtz _________________________________


Helmholtz :: Helmholtz (Mesh* m, boolean kioP) // kioP=false by default 
    : SemiDiscreteProblem(m) 
{
    // Constructor. Creates a problem on 'm'.

    alpha = UNINITIALIZED_REAL ;
    kioPressure = kioP;
}


Helmholtz :: ~Helmholtz ()
{
    // Destructor.
} 


real Helmholtz :: GetAlpha ()
{
    // Returns the coefficient alpha of the receiver.

    return alpha ;
}


void Helmholtz :: GetAx (Field* x, Field* y) 
{
    // Initializes 'y' to the operator of the receiver evaluated at 'x'.
    // The Dirichlet conditions are made homogeneous, consistently with the
    // formulation of the left-hand side (see method GetRHS).

    // set to 0 the columns of Qt.A.Q that are affected by bc's

    ApplyHomogeneousDirichletConditions(x) ;

    // x -> Q x (Q = mortar operator)
    x->Distribute() ;

    // y = A.Q x
    GetFlatAx(x->GetMain(),y->GetMain()) ;

    // y = Qt.A.Q x, assembled
    y->Assemble() ;

    // set to 0 the rows of Qt.A.Q that are affected by bc's
    ApplyHomogeneousDirichletConditions(y) ;
}

void Helmholtz :: GetFlatAx (FlatField* x, FlatField* y)
{
    // - if 'mode' = IMPLICIT, sets 'y' to the weak form of:
    //     lambda x - alpha Laplacian(x) ,
    //   i.e., to:
    //     lambda (x, psi) + alpha (grad x, grad psi)
    //
    // - if 'mode' = EXPLICIT, sets 'y' to the weak form of:
    //     lambda x ,
    //   i.e., to:
    //     lambda (x, psi)

# ifdef REQUIRE
    Require("same mesh", x->GetMesh() == mesh && y->GetMesh() == mesh) ;
    Require("same interpolation", y->HasSameInterpolationAs(unknown->GetMain()));
    Require("has integration rule", integrationRule != NULL) ;
    InMethod("Helmholtz::GetFlatAx(x,y)") ;
# endif

    FlatField *work ;

    if (alpha == 0 || mode == EXPLICIT) {
	if (lambda == 0)                              // y = 0
	    y->SetValues(ZERO) ;
	else {                                        // y = lambda B x
	    y->SetToWeakIdentity(x,integrationRule) ;
	    y->Multiply(lambda) ;
	}
    }

    else {
	y->SetToWeakLaplacian(x,integrationRule) ;
	y->Multiply(-alpha) ;                         // y = - alpha L x
	if (lambda != 0) {                            // y = (lambda B - alpha L) x
	    work = y->GetWork() ;
	    work->SetToWeakIdentity(x,integrationRule) ;
	    y->Add(work,lambda) ;
	    y->Retrieve(work) ;
	}
    }
}


void Helmholtz :: GetFlatRHS (FlatField* x)
{
    // Initializes 'x' to:
    //     (f,psi) + (g,psi) + (Hist_s, psi) + (Hist_w),
    // i.e., to
    //     B f + B_contour g + B Hist_s + Hist_w,
    // where Hist_s and Hist_w stand for the strong and the weak parts of the
    // history field, and B_contour is the mass operator on the contour mesh.
    //
    // Ignores the B f and B_contour g contributions if in explicit mode.

# ifdef REQUIRE
    Require("same mesh", x->GetMesh() == mesh) ;
    Require("has integration rule", integrationRule != NULL) ;
    InMethod("Helmholtz::GetFlatRHS(x)") ;
# endif

    FlatField *work ;
    int       i ;

    // 1. Body force f

    if (forcing && mode == IMPLICIT) {
	EvaluateForcing(time) ;
	x->SetToWeakIdentity(forcing,integrationRule) ;
    }
    else
	x->SetValues(ZERO) ;

    // 2. Contour force g

    if (mode == IMPLICIT)
	for (i=1 ; i <= neumannConditions->GetSize() ; i++)
	    neumannConditions->At(i)->ApplyOn(x,time) ;

    // 3. History vectors

    if (strongHistory) {
	work = unknown->GetMain()->GetWork() ;
	work->SetToWeakIdentity(strongHistory,integrationRule) ;
	x->Add(work) ;
	unknown->GetMain()->Retrieve(work) ;
    }

    if (weakHistory)
	x->Add(weakHistory) ;   
  
    if (kioPressure)
	x->SetToZeroAlgebraicMean(1) ;

    ////  printf(" END FLAT RHS: Neuclide %24.12e \n", x->GetEuclidianNorm());
}


void Helmholtz :: GetFunctional (FlatField* x, real t, FlatField* y)
{
    // Initializes 'y' to the weak form of: 
    //    alpha Laplacian(x) + f(t),
    // i.e, to:
    //   -alpha (grad x, grad psi) + (f(t), psi).
    // 't' may differ from the current time 'time'.

# ifdef REQUIRE
    Require("same mesh", x->GetMesh() == mesh && y->GetMesh() == mesh) ;
    Require("same interpolation", y->HasSameInterpolationAs(unknown->GetMain()));
    Require("has integration rule", integrationRule != NULL) ;
    InMethod("Helmholtz::GetFunctional(x,t,y)") ;
# endif

    FlatField *work ;

    y->SetToWeakLaplacian(x,integrationRule) ;
    y->Multiply(alpha) ;

    if (forcing) {
   
	EvaluateForcing(t) ;
	work = unknown->GetMain()->GetWork() ;
	work->SetToWeakIdentity(forcing,integrationRule) ;
	y->Add(work) ;
	unknown->GetMain()->Retrieve(work) ;
    }
}


void Helmholtz :: GetRHS (Field* rhs)
{
    // Initializes 'rhs' to the left-hand side of the operator of the receiver.
    // If 'rhs' is a mortared field, only its main field is concerned.
    //
    // The contributions of the non-homogeneous Dirichlet conditions are
    // accounted for here. The original system is:
    //           A u  = r
    // with DC:    u  = u0  on some of the values of u.
    // The linear system is modified, by splitting u into:
    //              u = u* + u0 ,
    // and thus reads (using the linearity of A):
    //           A u* = r - A u0 = r*
    //
    // This method initializes 'rhs' to r*, in assembled form.
    // The system can be solved as A u* = r*, with HOMOGENEOUS conditions on u*.
    // Then the solution of the original system is u = u* + u0.

    Field     *u0 ;
    FlatField *q ;

    GetFlatRHS(rhs->GetMain()) ;

    ////  printf(" RHS: b4 dirichletet assem  Neuclide %24.12e \n", rhs->GetEuclidianNorm());
    // get u0. if u (and thus u0) is a mortared field, the influence of the 
    // Dirichlet conditions on the interface values is calculated

    if (dirichletConditions) {

	u0 = unknown->FieldGetWork() ;
	u0->SetValues(ZERO) ;
	ApplyDirichletConditions(u0) ;
	u0->Distribute() ;
    
	// get r*

	q = unknown->GetMain()->GetWork() ;
	GetFlatAx(u0->GetMain(),q) ;

	////    printf(" q: Neuclide %24.12e \n", q->GetEuclidianNorm());
	rhs->GetMain()->Subtract(q) ;
    }
    // transfer the residual on true degrees of freedom, i.e., apply Q(transp)
    // to r* and assemble it
  
    ////  printf(" RHS: b4 assem Neuclide %24.12e \n", rhs->GetEuclidianNorm());

    rhs->Assemble() ;

    ////  printf(" RHS: AFT-ASSEM et Neumann Neuclide %24.12e \n", rhs->GetEuclidianNorm());
  
    // set the residual to 0 were Dirichlet bc's apply, in order to ensure that
    // the solution u* is 0 (homogeneous) there
  
    ApplyHomogeneousDirichletConditions(rhs) ;
  
    // housekeeping
    if (dirichletConditions) {
	unknown->Retrieve(u0) ;
	unknown->GetMain()->Retrieve(q) ;
    }
    ////  printf("FING  RHS: Neuclide %24.12e\n", rhs->GetEuclidianNorm());
}


boolean Helmholtz :: IsHelmholtz () 
{
    // Return 'true', because the receiver's type is Helmholtz.

    return true ;
}

void Helmholtz :: Print()
{
    printf(" Helmholtz :: Print():: alpha = %f  SDP: lambda = %f time = %f  ", alpha, lambda, time);
}

void Helmholtz :: SetAlpha (real val)
{
    // Initializes the alpha coefficient of the receiver to 'val'.

    alpha = val ;
}


void Helmholtz :: SubtractWeakHistory (FlatField* hist)
{
    // Adds "minus" 'hist' to the weak part (the one that does not need to be 
    // integrated) of the history field of the receiver.
    //
    // Remark. In order to speed up calculations and to save memory space, no 
    //         copy of 'hist' is made, so do not modify 'hist' before the 
    //         receiver is solved.

# ifdef REQUIRE
    Require("same mesh", hist->GetMesh() == mesh) ;
    Require("same interpolation", 
	    hist->HasSameInterpolationAs(unknown->GetMain())) ;
    InMethod("Helmholtz::SubtractWeakHistory(hist)") ;
# endif

    if (weakHistory)
	weakHistory->Subtract(hist) ;
    else {
	hist->SwitchSigns() ;
	ref(hist) ;
	weakHistory = hist ;
    }
}


//________________________________ SteadyStokes _______________________________


SteadyStokes :: SteadyStokes (Mesh* m, DecompositionType typ) 
    : SemiDiscreteProblem(m) 
{
    // Constructor. Creates a steady Stokes problem of type 'typ' defined on 'm'.

# ifdef REQUIRE
    Require("valid decomposition", typ == UZAWA || typ == BP1 || typ == BP3 ||
	    typ == BP1_PC || typ == BP3_PC || typ==YOSIDA ) ;
    InMethod("SteadyStokes::SteadyStokes(m,typ)") ;
# endif

    type                    = typ ;
    velocity                = NULL ;
    previousPressure        = NULL ;
    nu                      = UNINITIALIZED_REAL ;
    integrationRuleMomentum = NULL ;
    helmholtz               = NULL ;
    internalHelmholtz       = NULL ;
    vTilde                  = NULL ;

    stressfreeBC            = false;

    CreateHelmholtz () ;

    if (type == UZAWA)
	CreateInternalHelmholtz() ;
}


SteadyStokes :: ~SteadyStokes ()
{
    // Destructor.

    unref(helmholtz) ;
    unref(internalHelmholtz) ;
    unref(vTilde) ;
    unref(velocity) ;
    unref(previousPressure) ;
    unref(integrationRuleMomentum) ;
}


void SteadyStokes :: AddDirichletCondition (DirichletCondition*)
{
    // Issues an error message because pressure Dirichlet conditions are not
    // implemented.

    Error("SteadyStokes::AddDirichletCondition(dc)",
	  "Dirichlet conditions on the pressure are not allowed") ;
}


void SteadyStokes :: AddDirichletConditionV (DirichletCondition* dc)
{
    // Adds 'dc' to the list of velocity Dirichlet conditions of the receiver.

    DirichletCondition *dcHomog ;

//  dirichletConditionsV->Put(dc) ;

    helmholtz->AddDirichletCondition(dc) ;

    if (type == UZAWA) {
	dcHomog = dc->DuplicateHomogeneous() ;
	internalHelmholtz->AddDirichletCondition(dcHomog) ;
	unref(dcHomog) ;
    }
}


void SteadyStokes :: AddNeumannCondition (NeumannCondition*)
{
    // Issues an error message because pressure Neumann conditions are not
    // implemented.

    Error("SteadyStokes::AddNeumannCondition(nc)",
	  "Neumann conditions on the pressure are not allowed") ;
}


void SteadyStokes :: AddNeumannConditionV (NeumannCondition* nc)
{
    // Adds 'nc' to the list of velocity Neumann conditions of the receiver.

//  neumannConditionsV->Put(nc) ;

    helmholtz->AddNeumannCondition(nc) ;
}


void SteadyStokes :: ApplyHomogeneousDirichletConditionsV (Field* x)
{
    // Applies the velocity Dirichlet conditions of the receiver to 'x'. This is
    // done in a homogeneous way (the prescribed values are replaced by 0).

    helmholtz->ApplyHomogeneousDirichletConditions(x) ;
}


void SteadyStokes :: ApplyJ (MortaredField* x, MortaredField* y)
{
    // Applies y = J x, with homogeneous Dirichlet conditions.
    // In the Uzawa case, this amounts to solving H y = x. 
    // In the Blair Perot cases, this reduces to evaluating functions of x.
    // In transient analysis (unsteady Stokes) with explicit mode, Uzawa, BP and
    // BP_PC are identical.

# ifdef REQUIRE
    Require("if explicit mode, lambda != 0", mode==IMPLICIT || lambda != 0.) ;
    Require("has momentum integration rule", integrationRuleMomentum != NULL) ;
    InMethod("SteadyStokes:ApplyJ(x,y)") ;
# endif

    MortaredField *y1, *y2, *y3 ;

    if (type == UZAWA) {
	// J = H-1
	internalHelmholtz->SetUnknown(y) ;
	internalHelmholtz->SetWeakHistory(x->GetMain()) ; // use weakHistory, since
	internalHelmholtz->Solve() ;                      // 'x' is already in weak
	internalHelmholtz->SetUnknown(NULL) ;             // form
    }

    else if (type == BP1 || type == BP1_PC || type == YOSIDA ) {
	// J = (lambda B)-1
	y->GetMain()->CopyFrom(x->GetMain()) ;
	y->GetMain()->Multiply(ONE/lambda) ;
	helmholtz->ApplyInverseB(y) ;                   // use the B matrix defined
    }

    else if (type == BP3 || type == BP3_PC) {
	// J = c B-1 + c^2 nu (B-1 A) B-1 + c^3 nu^2 (B-1 A)^2 B-1
	// where c = 1/lambda.

	// y1 = B-1 x
	y1 = y ;
	y1->GetMain()->CopyFrom(x->GetMain()) ;
	helmholtz->ApplyInverseB(y1) ;

	if (mode == IMPLICIT) {
	    // y2 = (B-1 A) B-1 x
	    y2   = y->GetWork() ;
	    y2->GetMain()->SetToWeakLaplacian(y1->GetMain(),integrationRuleMomentum);
	    helmholtz->ApplyInverseB(y2) ;

	    // y3 = (B-1 A)^2 B-1 x
	    y3 = y->GetWork() ;
	    y3->GetMain()->SetToWeakLaplacian(y2->GetMain(),integrationRuleMomentum) ;
	    helmholtz->ApplyInverseB(y3) ;
	}

	// y = c y1 - c^2 nu y2 + c^3 nu^2 y3
	y->GetMain()->Multiply(ONE/lambda) ;
	if (mode == IMPLICIT) {
	    y->GetMain()->Add(y2->GetMain(),nu/(lambda*lambda)) ;
	    y->GetMain()->Add(y3->GetMain(),nu*nu/(lambda*lambda*lambda)) ;
	    y->Retrieve(y2) ;
	    y->Retrieve(y3) ;
	}
    }

    else
	Error("SteadyStokes::ApplyJ(x,y)","unexpected decomposition type") ;
}


void SteadyStokes :: CreateHelmholtz ()
{
    // Creates the Helmholtz problem for solving the momentum equation (i.e., for
    // computing the tentative velocity v~ in phase 1).

    Solver *s ;

    helmholtz = new Helmholtz(mesh) ;
    helmholtz->SetPreconditioner(NULL) ;
    helmholtz->SetTime(time) ;

    // default solver (usually overridden by the user).
    s = new ConjugateGradient(500,DEFAULT_TOLERANCE) ;
    s->SetName("Helmholtz") ;
    s->SetVerbosity(0) ; 
    helmholtz->SetSolver(s) ; 
    unref(s) ;
}


void SteadyStokes :: CreateInternalHelmholtz ()
{
    // Creates the Helmholtz problem for solving the Helmholtz equations in
    // phases 2 (pressure) and 3 (correction of the tentative velocity).
 
    Solver *s ;

    internalHelmholtz = new Helmholtz(mesh) ;
    internalHelmholtz->SetPreconditioner(NULL) ;

    // default solver (usually overridden by the user). Note that a higher 
    // precision is required in the internal loop than in the main pressure loop
    // of phase 2

    s = new ConjugateGradient(200,DEFAULT_TOLERANCE/100.) ;
    s->SetName("internal") ;
    s->SetVerbosity(0) ;
    internalHelmholtz->SetSolver(s) ;
    unref(s) ;
}
FullMatrix* SteadyStokes :: GetMatrixOperatorV(int q, real coeff )
{
    // Returns the matrix operator A of the receiver. Constructs it by 'probing'
    // (expensive!).
    // The order of A is equal to the number of degrees of freedom of the mesh.
    // Boundary conditions and mortar constraints are taken into account.

    MortaredField    *temp1, *temp2, *work , *work2 ;
    FullMatrix *A ;
    int        i, nbdof ;

    // 1. Compute A
    // 1.a Initialisation 

    printf(" Enter FullMatrix* SteadyStokes :: GetMatrixOperatorV(int q, real coeff ) \n");
    nbdof = velocity->GetNbDofs() ;

    temp1 = velocity->GetWork();
    temp2 = velocity->GetWork(1);
    work  = velocity->GetWork(1);
    work2  = velocity->GetWork();

    printf(" Nbdofs = %d \n",nbdof);

    A     = new FullMatrix(nbdof,nbdof) ;

    for (i=1 ; i<=nbdof ; i++)
	for (int j=1 ; j<=nbdof ; j++)    
	    A->At(i,j) = ZERO;

    for (i=1 ; i<=nbdof ; i++) {
	if (i>1) 
	    temp1->SetDof(i-1,ZERO) ;
	temp1->SetDof(i,ONE) ;
	ApplyHomogeneousDirichletConditions(temp1) ;
	temp1->Distribute() ;

	// 1.b Set temp1 as "previous" value in time schemes

	temp1->GetMain()->Multiply(coeff);

	helmholtz->SetStrongHistory(temp1->GetMain());
	helmholtz->Solve() ;

	work->GetMain()->SetToWeakDivergence(velocity->GetMain(),integrationRule);
	GetAx(work,temp2);
	temp2->Assemble() ;
//     ApplyHomogeneousDirichletConditions(temp2) ;
	temp2->Distribute();

	work2->GetMain()->SetToWeakGradient(temp2->GetMain(),integrationRule);

//     ApplyHomogeneousDirichletConditions(work2) ;
//     work2->Distribute() ;
	//    work2->GetMain()->SwitchSigns();

	helmholtz->SetUnknown(velocity);
    
	temp1->Multiply(-1.);
	helmholtz->SetStrongHistory(temp1->GetMain());
	helmholtz->SetWeakHistory(work2->GetMain());
	helmholtz->Solve() ;

	velocity->Assemble() ;
	ApplyHomogeneousDirichletConditions(velocity) ;
	velocity->Distribute();

	velocity->CopyDofsToMatrixColumn(A,i) ;

    }
  
    // 3. Take homogeneous Dirichlet bc's into account in LHS, using a trick for
    //    detecting where Dirichlet bc's apply

    temp1->SetValues(ONE) ;
    ApplyHomogeneousDirichletConditions(temp1) ;
    temp1->Distribute() ;

    for (i=1 ; i<=nbdof ; i++)
	if (temp1->GetDof(i) != ONE)
	    A->At(i,i) = ONE ;

    velocity->Retrieve(temp1) ;
    velocity->Retrieve(temp2) ;
    velocity->Retrieve(work) ;
    velocity->Retrieve(work2) ;

    return A ;
}


void SteadyStokes :: GetAx (Field* x, Field* y)
{
    // Initializes 'y' to 'A x', where A is the LHS of phase 2, i.e., the LHS of:
    //   -D J D(T) x = D v~ .

# ifdef REQUIRE
    Require("x is a flat field", x->IsFlatField()) ;
    Require("y is a flat field", y->IsFlatField()) ;
    Require("has mass integration rule", integrationRule != NULL) ;
    InMethod("SteadyStokes::GetAx(x,y)") ;
# endif

    MortaredField *work1, *work2 ;
    FlatField     *xx, *yy ;

    xx   = x->GetMain() ;
    yy   = y->GetMain() ;

    // seems to bug
    // xx->SetToZeroAlgebraicMean(1) ;   // squeeze out any component of x in the 
    // null space of D J D(T)
    // 1. D(T) x

    work1 = velocity->GetWork() ;
    work1->GetMain()->SetToWeakGradient(xx,integrationRule) ;

    // 2. J D(T) x

    work2 = velocity->GetWork() ;
    ApplyJ(work1,work2) ;
    velocity->Retrieve(work1) ;
 
    // 3. D J D(T) x

    yy->SetToWeakDivergence(work2->GetMain(),integrationRule) ;
    velocity->Retrieve(work2) ;

}


real SteadyStokes :: GetC1 ()
{
    // Returns the c1 coefficient of the Couzy preconditioner
    //   inv(P) = c1 inv(B) + c2 inv(E).

    if (type == UZAWA)
	return nu ;
    else                          // BP or BP_PC
	return ZERO ;
}


real SteadyStokes :: GetC2 ()
{
    // Returns the c2 coefficient of the Couzy preconditioner
    //   inv(P) = c1 inv(B) + c2 inv(E).
    // Remark. Using a Couzy preconditioner in the Uzawa case is uneffective if
    //         lambda = 0 (steady case); using a mass preconditioner is 
    //         equivalent and cheaper.

    return lambda ;
}


void SteadyStokes :: GetFlatAx (FlatField* x, FlatField* y)
{
    // Same as method GetAx.
    // Typically used for computing the diagonal of the pressure operator.

    GetAx(x,y) ;
}


void SteadyStokes :: GetFunctional (FlatField* x, real t, FlatField* y)
{
    // Initializes 'y' to the weak form of: 
    //    - grad p + nu Laplacian(x) + f(t),
    // The current value of the pressure is used for 'p'.

# ifdef REQUIRE
    Require("same mesh", x->GetMesh() == mesh && y->GetMesh() == mesh) ;
    Require("same interpolation",y->HasSameInterpolationAs(velocity->GetMain()));
    Require("has momentum integration rule", integrationRuleMomentum != NULL) ;
    InMethod("SteadyStokes::GetFunctional(x,t,y)") ;
# endif

    FlatField *work, *f ;

    work = velocity->GetMain()->GetWork() ;

    y->SetToWeakLaplacian(x,integrationRuleMomentum) ;
    y->Multiply(nu) ;

    work->SetToWeakGradient(unknown->GetMain(),GetIntegrationRuleMass()) ;
    y->Subtract(work) ;

    f = helmholtz->GetForcing() ;
    if (f) {
	helmholtz->EvaluateForcing(t) ;
	work->SetToWeakIdentity(f,integrationRuleMomentum) ;
	y->Add(work) ;
    }

    velocity->GetMain()->Retrieve(work) ;
}


Helmholtz* SteadyStokes :: GetHelmholtz ()
{
    // A debugging tool.
    // Returns the Helmholtz solver (for the system H v~ = B f + B_contour g) of
    // the receiver.

    return helmholtz ;
}


Vector<DirichletCondition*>* SteadyStokes :: 
GetHomogeneousDirichletConditionsV ()
{
    // Returns a new vector containing copies of the velocity Dirichlet 
    // conditions of the receiver, with all coefficients of the prescribed-value 
    // fields set to 0.

    Vector<DirichletCondition*> *dcs, *answer ;
    int                         i ;

    answer = new Vector<DirichletCondition*>() ;

    dcs = helmholtz->GetDirichletConditions() ;
    for (i=1 ; i <= dcs->GetSize() ; i++)
	answer->Put(dcs->At(i)->DuplicateHomogeneous()) ;

    return answer ;
}


FlatField* SteadyStokes :: GetIntegrationRuleMass ()
{
    // Returns the integration rule for the mass equation.
    // Unless specified otherwise, it is identical to the pressure field.

    return integrationRule ;
}


FlatField* SteadyStokes :: GetIntegrationRuleMomentum ()
{
    // Returns the integration rule for the momentum equation.
    // Unless specified otherwise, it is identical to the velocity field.

    return integrationRuleMomentum ;
}


Helmholtz* SteadyStokes :: GetInternalHelmholtz ()
{
    // A debugging tool.

    return internalHelmholtz ;
}


real SteadyStokes :: GetNu ()
{
    // Returns the viscosity of the receiver.

    return nu ;
}

FlatField*  SteadyStokes :: GetPreviousPressure()
{
    return previousPressure;
}

void SteadyStokes :: GetRHS (Field* rhs)
{
    // Initializes 'rhs' to the RHS of phase 2, i.e., the RHS of: 
    //   -D J D(T) p = D v~.

# ifdef REQUIRE
    Require("'rhs' is a flat field", rhs->IsFlatField()) ;
    Require("'rhs' has 1 component", rhs->GetNbComponents() == 1) ;
    Require("has mass integration rule", integrationRule != NULL) ;
    InMethod("SteadyStokes::GetRHS(rhs)") ;
# endif

    ((FlatField*)rhs)->SetToWeakDivergence(velocity->GetMain(),integrationRule) ;

    // make sure that rhs has no component in the null space of Uzawa operator
    // (even a tiny such component may prevent the conjugate-gradient iterations
    // for the pressure from converging)

    if (!stressfreeBC)
	rhs->GetMain()->SetToZeroAlgebraicMean(1) ;

}


boolean SteadyStokes :: IsSteadyStokes () 
{
    // Return 'true', because the receiver's type is SteadyStokes.

    return true ;
}


void SteadyStokes :: SetDirichletConditions (Vector<DirichletCondition*>*)
{
    // Issues an error message because pressure Dirichlet conditions are not
    // implemented.

    Error("SteadyStokes::SetDirichletConditions(dc)",
	  "Dirichlet conditions on the pressure are not allowed") ;
}


void SteadyStokes :: SetDirichletConditionsV (Vector<DirichletCondition*>* dcs)
{
    // Sets the set of velocity Dirichlet conditions of the receiver to 'dcs'.

    Vector<DirichletCondition*> *homog ;
    int                         i ;

    // the Helmholtz problem uses the bc's
    helmholtz->SetDirichletConditions(dcs) ;

    // the internal Helmholtz problem uses the homogeneous bc's
    if (type == UZAWA) {
	homog = GetHomogeneousDirichletConditionsV() ;
	internalHelmholtz->SetDirichletConditions(homog) ;
	for (i=1 ; i <= homog->GetSize() ; i++)
	    unref(homog->At(i)) ;
	delete homog ;
    }
}


void SteadyStokes :: SetForcing (FlatField* f)
{
    // Sets the forcing term of the receiver to 'f'.
    // 'f' may be NULL.

# ifdef REQUIRE
    Require("same mesh", f == NULL || f->GetMesh() == mesh) ;
    InMethod("SteadyStokes::SetForcing(f)") ;
# endif

    helmholtz->SetForcing(f) ;
}


void SteadyStokes :: SetForcingFunction (Function* funct)
{
    // Sets the time function for the forcing term of the receiver to 'funct'.
    // 'funct' may be NULL.

    helmholtz->SetForcingFunction(funct) ;
}


void SteadyStokes :: SetIntegrationRuleMass (FlatField* rule)
{
    // Initializes the integration rule for the mass equation to 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("SteadyStokes::SetIntegrationRuleMass(rule)") ;
# endif

    SetIntegrationRule(rule) ;
}


void SteadyStokes :: SetIntegrationRuleMomentum (FlatField* rule)
{
    // Initializes the integration rule for the momentum equations to 'rule'.

# ifdef REQUIRE
    Require("same mesh", rule->GetMesh() == mesh) ;
    Require("'rule' has 1 component", rule->GetNbComponents() == 1) ;
    InMethod("SteadyStokes::SetIntegrationRuleMomentum(rule)") ;
# endif

    ref(rule) ;
    unref(integrationRuleMomentum) ;

    integrationRuleMomentum = rule ;

    helmholtz->SetIntegrationRule(rule) ;

    if (type == UZAWA)
	internalHelmholtz->SetIntegrationRule(rule) ;
}


void SteadyStokes :: SetLambda (real val)
{
    // Sets the lambda coefficient of the receiver to 'val'.

# ifdef REQUIRE
    Require("lambda non nul if BP decomposition", type == UZAWA || val != 0.) ;
    InMethod("SteadyStokes::SetLambda(val)") ;
# endif

    lambda = val ;

    helmholtz->SetLambda(val) ;

    if (type == UZAWA)
	internalHelmholtz->SetLambda(val) ;
}


void SteadyStokes :: SetNeumannConditions (Vector<NeumannCondition*>*)
{
    // Issues an error message because pressure Neumann conditions are not
    // implemented.

    Error("SteadyStokes::SetNeumannConditions(nc)",
	  "Neumann conditions on the pressure are not allowed") ;
}


void SteadyStokes :: SetNeumannConditionsV (Vector<NeumannCondition*>* ncs)
{
    // Sets the set of velocity Neumann conditions of the receiver to 'ncs'.

    helmholtz->SetNeumannConditions(ncs) ;
}

void SteadyStokes :: ImposeNoPressureNullSpace()
{
    stressfreeBC = true;
}


void SteadyStokes :: SetNu (real val)
{
    // Sets the kinematic viscosity (nu) of the receiver to 'val'.
  
    nu = val ;

    helmholtz->SetAlpha(val) ;

    if (type == UZAWA)
	internalHelmholtz->SetAlpha(val) ;
}


void SteadyStokes :: SetPressure (FlatField* p)
{
    // Sets the pressure field of the receiver to 'p'.
    // In the case of a BP_PC decomposition, 'p' provides the values of p(0).

# ifdef REQUIRE
    Require("same mesh", p->GetMesh() == mesh) ;
    InMethod("SteadyStokes::SetPressure(p)") ;
# endif

    SetUnknown(p) ;

    if (type == BP1_PC || type == BP3_PC || type==YOSIDA) {
	unref(previousPressure) ;
	previousPressure = p->DuplicateEmpty() ;
    }
}


void SteadyStokes :: SetStrongHistory (FlatField* x)
{
    // Sets the "strong" part of the history field of the receiver to 'x'.
    // Applies to the time-differential part of the problem, i.e., to phase 1.

    helmholtz->SetStrongHistory(x) ;
} 


void SteadyStokes :: SetTemporalUnknown (Field* v)
{
    // Sets the velocity field of the receiver to 'v'.
    // Remark. In this class the temporal unknown and the unknown are not the
    //         same field (velocity vs pressure).

# ifdef REQUIRE
    Require("same mesh", v->GetMesh() == mesh) ;
    Require("'v' is a mortared field", v->IsMortaredField()) ;
    InMethod("SteadyStokes::SetTemporalUnknown(v)") ;
# endif

    SetVelocity((MortaredField*)v) ;
}


void SteadyStokes :: SetTime (real t)
{
    // Sets the clock of the receiver to 't'. Used when evaluating the forcing
    // term and the time-dependent Dirichlet conditions.

    time = t ;

    helmholtz->SetTime(t) ;
}


void SteadyStokes :: SetToExplicit ()
{
    // Whenever the LHS and RHS operators of the momentum equation are 
    // calculated, ignores all terms of the functional F(v).

    mode = EXPLICIT ;

    helmholtz->SetToExplicit() ;

    if (type == UZAWA)
	internalHelmholtz->SetToExplicit() ;
}


void SteadyStokes :: SetToImplicit ()
{
    // Whenever the LHS and RHS operators of the momentum equation are 
    // calculated, considers all terms of the functional F(v).
    // This is the default mode.

    mode = IMPLICIT ;

    helmholtz->SetToImplicit() ;

    if (type == UZAWA)
	internalHelmholtz->SetToImplicit() ;
}


void SteadyStokes :: SetVelocity (MortaredField* v)
{
    // Initializes the velocity field of the receiver to 'v'.

# ifdef REQUIRE
    Require("same mesh", v->GetMesh() == mesh) ;
    Require("same interpolation", velocity == NULL || 
	    v->HasSameInterpolationAs(velocity)) ;
    InMethod("SteadyStokes::SetVelocity(v)") ;
# endif

    ref(v) ;
    unref(velocity) ;
    velocity = v ;

    helmholtz->SetUnknown(v) ;

    if (type == UZAWA)
	internalHelmholtz->SetUnknown(v) ;    // used only for interpolation
}


void SteadyStokes :: SetWeakHistory (FlatField* x)
{
    // Sets the "weak" part of the history field of the receiver to 'x'.
    // Applies to the time-differential part of the problem, i.e., to phase 1.

    helmholtz->SetWeakHistory(x) ;
} 


void SteadyStokes :: Solve ()
{
    // Solves for the velocity and the pressure at the current time step, using
    // the generalized LU decomposition method.

# ifdef REQUIRE
    Require("has mass integration rule", integrationRule != NULL) ;
    InMethod("SteadyStokes::Solve()") ;
# endif

    MortaredField *work1, *work2 ;
    FlatField     *pressure, *Gp ;

    FlatField     *divVtildeOnPmesh, *divVtildeOnVmesh; // GODA

    // Phase 1: compute the tentative velocity v~ by: 
    //             H v~ = B f + B_contour g  (+ D(T) pN if BP_PC decoupling), 
    //          with nonhomog bc's. v~ overwrites v.

//cout_root << "[MONITORING] PHASE 1" << endl;
//	_mpimon_start_display();

    pressure = (FlatField*)unknown ;

    if (type == BP1_PC || type == BP3_PC || type == YOSIDA) {
	previousPressure->CopyFrom(pressure) ;
	Gp = velocity->GetMain()->GetWork() ;                        
	Gp->SetToWeakGradient(previousPressure,
			      GetIntegrationRuleMass()) ;      // -D(T) pN
	helmholtz->SubtractWeakHistory(Gp);
    }

    helmholtz->Solve() ;

#ifdef FILTER_VTILDE
    velocity->GetMain()->Filter(0.05, 1);
#endif

    if (type == BP1_PC || type == BP3_PC || type == YOSIDA)
	velocity->GetMain()->Retrieve(Gp) ;

// # ifdef CHECK
//   unref(vTilde) ;                             // useful to test if pressure
//   vTilde = velocity->GetMain()->Duplicate() ; // operator correctly implemented
//   printf("v~: "); velocity->GetMain()->Print();
// # endif


    // Phase 2: compute the pressure p: -D J D(T) p = D v~
    //          (if BP_PC, delta_p = p-pN is obtained, not p)
    //          uses the homogeneous bc's for J

//_mpimon_display_counters();
//cout_root << "[MONITORING] PHASE 2" << endl;

    if (internalHelmholtz)
	internalHelmholtz->GetSolver()->SetAdjustTolerance(true) ;


    GetSolver()->Solve(this) ;
    pressure = (FlatField*)unknown ;

    if (internalHelmholtz)
	internalHelmholtz->GetSolver()->SetAdjustTolerance(false) ;


    // Phase 3: correct the tentative velocity v~ into v = v~ + J D(T) p .
    //          uses the homogeneous bc's for J.
    //          note the sign: WeakGradient(p) = - D(T) p

//_mpimon_display_counters();
//cout_root << "[MONITORING] PHASE 3" << endl;

    if (type == YOSIDA) {
	Gp = velocity->GetMain()->GetWork() ;                        
	Gp->SetToWeakGradient(pressure, GetIntegrationRuleMass()) ;

	helmholtz->SubtractWeakHistory(Gp);
	helmholtz->Solve() ;
	velocity->GetMain()->Retrieve(Gp) ;

    }
    else
    {
	work1 = velocity->GetWork() ;
	work2 = velocity->GetWork() ;
      
	if (type != UZAWA ){// GODA
	    // computes the additional term to the pressure to get an overall
	    // consistent scheme (cf. Timmermans's thesis)
	    divVtildeOnPmesh = pressure->GetWork(1) ;
	    divVtildeOnVmesh = velocity->GetMain()->GetWork(1) ;
	    divVtildeOnVmesh->SetToDivergence(velocity->GetMain(),true);
	    divVtildeOnPmesh->CopyInterpolateFrom(divVtildeOnVmesh,1,1);
	    divVtildeOnPmesh->Multiply(ONE_MINUS*this->GetNu());
	}
      
	work1->GetMain()->SetToWeakGradient(pressure,GetIntegrationRuleMass()) ;
	ApplyJ(work1,work2) ;
	velocity->GetMain()->Subtract(work2->GetMain()) ;
      
	velocity->Retrieve(work1) ;
	velocity->Retrieve(work2) ;
    }
    // Complementary phase:
    //          If BP_PC, get p = pN + delta_p
    //          Set mean of p to 0.

//_mpimon_display_counters();
//cout_root << "[MONITORING] PHASE COMPLEMENTAIRE" << endl;
  
    if (type == BP1_PC || type == BP3_PC || type == YOSIDA)
	pressure->Add(previousPressure) ;

    if (type != UZAWA){ // GODA
	// Adds additional term to the pressure to get consistent scheme
	if (type == BP1_PC || type == BP3_PC )
	    pressure->Add(divVtildeOnPmesh);
	pressure->Retrieve(divVtildeOnPmesh) ;
	velocity->GetMain()->Retrieve(divVtildeOnVmesh) ;
    }

    if (!stressfreeBC)
	pressure->SetToZeroMean(1) ;

#ifdef FILTER
    velocity->GetMain()->Filter(0.05, 1);
#endif
//_mpimon_display_counters();
//_mpimon_stop_display();

}


