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

// File main1D_heat.cxx
//
// Solves the 1D heat equation:
//
//     u,  = alpha u,   + f(x,t)
//       t           xx
// with:
//   . x in ]0,20[ and t >= 0
//   . alpha = 200 / (pi*pi)
//   . f(x,t) = sin(pi x / 10) exp(-t)
//   . boundary conditions: u(0,t) = u(20,t) = 0.
//   . initial conditions:  u(x,0) = sin(pi x / 20).
//
// The exact solution is:
//     u(x,t) = sin (pi x /10) exp(-t)
// and is identical to f(x,t).
//
// The measured error is the H1 norm of u(x,t) - u_exact(x,t), t given, where
// u_exact is evaluated on the same mesh as u(x,t). This implies that this
// measure of the error tends to 0 when the time increment delta_t tends to 0,
// although the space discretization is not accurate.


// Define (resp. undefine) USE_CG if you want a conjugate-gradient solver
// (resp. a direct solver) for solving the linear system. The output is the
// same. The standard solver is the C.G. one; however, by using a direct 
// solver the user can figure out more easily how the solution phase works.
#define USE_CG

// Choose (i.e., define) ONE among the following time-integration schemes:
//#define FORWARD_EULER
//#define BACKWARD_EULER
#define CRANK_NICOLSON
//#define ADAMS_BASHFORTH_1           // identical to FORWARD_EULER
//#define ADAMS_BASHFORTH_2
//#define ADAMS_BASHFORTH_3
//#define BDF_1                       // identical to BACKWARD_EULER
//#define BDF_2
//#define BDF_3
//#define BDF_4
//#define RUNGE_KUTTA_2
//#define RUNGE_KUTTA_4


#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. number of elements
  //   2. polynomial degree of the GLL rule
  //   3. number of time steps (default=1)
  //   4. time increment delta_t (default=1.).

  real dt ;
  int  nelem, polydeg, nstep ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 3 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      nstep   = 1 ;
      dt      = 1. ;
      break ;
    case 4 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      nstep   = atoi(argv[3]) ;
      dt      = 1. ;
      break ;
    case 5 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      nstep   = atoi(argv[3]) ;
      dt      = atof(argv[4]) ;
      break ;
    default :
      printf("arguments are: nelem, polydeg, nstep, dt.\n\n") ;
      exit(0) ;
      break ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2 ;
  Vertex       *v1, *v2 ;
  Edge         *line ;
  Distribution *distrib ;
  Mesh1D       *mesh ;

  p1 = new Point( 0.) ;
  p2 = new Point(20.) ;

  v1 = new Vertex(p1) ;
  v2 = new Vertex(p2) ;

  line = new Line(v1,v2) ;

  distrib = new UniformDistribution(nelem) ;
  line->SetDistribution(1,distrib) ;

  mesh = line->GenerateMesh() ;
  mesh->DispatchElements() ;
  //  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent ;
  FlatField     *coord ;
  MortaredField *vel ;

  //    1. Parent element definitions
  //    -------------------------------

  parent  = ParentElementGLL::GetParentElementOfDegree(polydeg) ;
  mesh->SetDefaultParentElement(parent) ;

  //    2. Field declarations and initializations
  //    -------------------------------------------

  coord = mesh->GetCoordinates() ;
  //  printf("\nCoordinates :") ; coord->Print() ;

  vel = new MortaredField(1,coord) ;


  // ============================================
  // *** Forcing term and boundary conditions ***
  // ============================================

  Mesh               *boundary ;
  FlatField          *forcing, *prescribed ;
  Function           *forcingFunction, *initFunction ;
  DirichletCondition *dc ; 

  //    1. Forcing term
  //    ---------------

  forcing         = vel->GetMain()->DuplicateEmpty() ;
  forcingFunction = new SinExp(PI/10.,-1.) ;
  //  printf("\nForcing term (initial): ") ; forcing->Print() ;

  //    2. Boundary conditions
  //    ----------------------

  boundary = mesh->GetSkin() ;
  //  printf("\nBoundary: ") ; boundary->Print() ;

  prescribed = new FlatField(1,coord,boundary) ;
  //  printf("\nPrescribed: "); prescribed->Print() ;

  dc = new DirichletCondition(prescribed,1) ;

  //    3. Initial conditions
  //    ---------------------

  initFunction = new SinExp(PI/10, 0.) ;               // f(x) = sin (pi x /10)
  vel->SetTo(initFunction) ;
  //  printf("\nInitial values: ") ; vel->Print() ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *rule ;

  rule = vel->GetMain()->DuplicateEmpty(1) ;


  // ===============
  // *** Problem ***
  // ===============

  HeatEquation          *problem ;
  TimeIntegrationScheme *scheme ;
  Solver                *solver ;

  //    1. Problem
  //    ----------

  problem = new HeatEquation(mesh) ;
  problem->SetThermalDiffusivity(200./(PI*PI)) ;
  problem->SetField(vel) ;
  problem->AddDirichletCondition(dc) ;
  problem->SetIntegrationRule(rule) ;
  problem->SetForcing(forcing) ;
  problem->SetForcingFunction(forcingFunction) ;


  //    2. Time-integration scheme
  //    --------------------------

  scheme = NULL ;

# ifdef FORWARD_EULER
    scheme = new ThetaMethod(0.) ;
# endif

# ifdef BACKWARD_EULER
    scheme = new ThetaMethod(1.) ;
# endif

# ifdef CRANK_NICOLSON
    scheme = new ThetaMethod(0.5) ;
# endif

# ifdef ADAMS_BASHFORTH_1
    scheme = new AdamsBashforth(1) ;
# endif

# ifdef ADAMS_BASHFORTH_2
    scheme = new AdamsBashforth(2) ;
# endif

# ifdef ADAMS_BASHFORTH_3
    scheme = new AdamsBashforth(3) ;
# endif

# ifdef BDF_1
    scheme = new BDF(1);
# endif

# ifdef BDF_2
    scheme = new BDF(2);
# endif

# ifdef BDF_3
    scheme = new BDF(3);
# endif

# ifdef BDF_4
    scheme = new BDF(4);
# endif

# ifdef RUNGE_KUTTA_2
    scheme = new RungeKutta(2);
# endif

# ifdef RUNGE_KUTTA_4
    scheme = new RungeKutta(4);
# endif

  scheme->SetTimeIncrement(dt) ;
  problem->SetTimeIntegrationScheme(scheme) ;

  //    3. Solver
  //    ---------

# ifdef USE_CG
    solver = new ConjugateGradient(300,1e-12) ;
# else
    solver = new DirectSolver() ;
# endif

  problem->SetSolver(solver) ;


  // ===================================
  // *** Solution and postprocessing ***
  // ===================================

  FlatField *error ;
  Function  *exactFunction ;
  FILE      *file ;
  real      norm ;
  int       i ;

  error         = vel->GetMain()->DuplicateEmpty() ;
  exactFunction = new SinExp(PI/10.,-1.) ;

  for (i=1 ; i<=nstep ; i++) {

    // approximate solution v
    if (i%1 == 0) printf("\n\nStart solving step %d\n\n",i) ;
    problem->Solve(1) ;

    // exact solution v_exact
    error->SetTo(exactFunction,scheme->GetTime()) ;

    // error = v - v_exact
    error->Subtract(vel->GetMain()) ;
    norm = error->GetH1Norm() ;
    if (i%1000 == 0) 
      printf("\nH1 norm of error: % .4e\n\n",norm) ;
  }

  printf("\nH1 norm of error on solution of step %d: %.4e\n\n",i-1,norm) ;

  file = Fopen("heat.dat","a+") ;
  fprintf(file,"%e %e\n",dt,norm) ;
  fclose(file) ;


  // ================
  // *** Cleaning ***
  // ================

  unref(boundary) ;
  unref(dc) ;
  unref(distrib) ;
  unref(error) ;
  unref(exactFunction) ;
  unref(forcing) ;
  unref(forcingFunction) ;
  unref(initFunction) ;
  unref(line) ;
  unref(mesh) ;
  unref(p1) ;
  unref(p2) ;
  unref(prescribed) ;
  unref(problem) ;
  unref(rule) ;
  unref(solver) ;
  unref(scheme) ;
  unref(v1) ;
  unref(v2) ;
  unref(vel) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}
