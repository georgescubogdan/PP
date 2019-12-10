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

// main2D_navier.cxx
//
// Solves the following 2D unsteady Navier-Stokes problem:
//
//     dv/dt + v.grad v = - grad p + 1/Re Laplacian v + f(x,y,t)
//     div v            = 0
//
//     Exact solution: u(x,y,t) =   - cx sy st
//                     v(x,y,t) =     sx cy st
//                     p(x,y,t) = -pi sx sy st + constant
//     with the notations:
//           cx = cos (pi/2 x), sx = sin (pi/2 x),
//           cy = cos (pi/2 y), sy = sin (pi/2 y),
//           ct = cos (pi/2 t), st = sin (pi/2 t).
//
// The data are:
//   . space domain: the square cavity [-1,1] x [-1,1] modelled by 2x2 
//     elements;
//   . Re = 1;
//   . f(x,y,t) = a complex expression given in class ForcingNavierStokes2D;
//   . Dirichlet conditions (on the velocity only): at all times, u and v are
//     set to the exact values on the 4 edges.
//   . initial conditions: at t=0, v=0 on the whole cavity (= exact values);
//     idem for the values at t=-dt, t=-2*dt, for the initial steps of 2-step,
//     3-step and 4-step time-integration schemes. 
//
// Karniadakis's method BDFn/EXn (n=1, 2, 3 or 4) is used. See Speculoos
// programmer's manual or Couzy's thesis, p 40.


// In order to choose the time-integration scheme to be used for the diffusive
// part of the problem, define (i.e., uncomment) ONE among the following
// options:

//#define BDF1
//#define BDF2
//#define BDF3
#define BDF4


// Define (i.e., uncomment) ONE among the following lines (the most accurate
// is the Uzawa one):

#define DECOMPOSITION_TYPE UZAWA
//#define DECOMPOSITION_TYPE BP1
//#define DECOMPOSITION_TYPE BP3
//#define DECOMPOSITION_TYPE BP1_PC
//#define DECOMPOSITION_TYPE BP3_PC


// Define the following option if you want the solution fields to be written
// into a file for Tecplot visualization:

//#define POSTPROC



#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //   2. nelem2:   number of elements in the Y direction.
  //   3. polydeg1: polynomial degree of the GLL rule for the velocity field in
  //                the X direction, in all elements.
  //   4. polydeg2: polynomial degree of the GLL rule for the velocity field in
  //                the Y direction, in all elements.
  //   5. nsteps:   the number of time steps;
  //   6. dt:       time increment delta_t;

  real dt ;
  int  nelem1, nelem2, polydeg1, polydeg2, nsteps ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem1   = 4 ; 
      nelem2   = 4 ; 
      polydeg1 = 7 ;
      polydeg2 = 7 ;
      nsteps   = 30 ;
      dt       = 0.033333333333333333333333333333 ;
      break ;
    case 7 :
      nelem1   = atoi(argv[1]) ;
      nelem2   = atoi(argv[2]) ;
      polydeg1 = atoi(argv[3]) ;
      polydeg2 = atoi(argv[4]) ;
      nsteps   = atoi(argv[5]) ;
      dt       = atof(argv[6]) ;
      break ;
    default :
      printf("%s%s%s", "valid number of arguments:\n",
             "  0 or\n",
             "  6 (nelem1,nelem2,polydeg1,polydeg2,nsteps,dt).\n\n") ;
      exit(0) ;
      break ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2, *p3, *p4 ;
  Vertex       *v1, *v2, *v3, *v4 ;
  Edge         *edge1, *edge2, *edge3, *edge4 ;
  Face         *face ;
  Distribution *distrib1, *distrib2 ;
  Mesh2D       *mesh ;

  p1 = new Point(0.,0.) ;
  p2 = new Point(1.,0.) ;
  p3 = new Point(1.,1.) ;
  p4 = new Point(0.,1.) ;

  v1 = new Vertex(p1) ;
  v2 = new Vertex(p2) ;
  v3 = new Vertex(p3) ;
  v4 = new Vertex(p4) ;

  edge1 = new Line(v1,v2) ;
  edge2 = new Line(v2,v3) ;
  edge3 = new Line(v3,v4) ;
  edge4 = new Line(v4,v1) ; 

  face  = new Quad(edge1,edge3,edge4,edge2) ;

  nelem1   = 2 ;
  nelem2   = 2 ;
  distrib1 = new UniformDistribution(nelem1) ;
  distrib2 = new UniformDistribution(nelem2) ;
  face->SetDistribution(1,distrib1) ;
  face->SetDistribution(2,distrib2) ;

  mesh = face->GenerateMesh() ;
  mesh->DispatchElementsByBlocks(nelem1,nelem2) ;
  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement          *parent1, *parent2 ;
  FlatField              *coord, *pressure ;
  MortaredField          *velocity ;
  Vector<ParentElement*> *parents ;
  int                    polydegP, shift ;

  //    1. Parent element definitions and coordinate field
  //    --------------------------------------------------

  parent1 = ParentElementGLL::GetParentElementOfDegree(polydeg1) ;
  parent2 = ParentElementGLL::GetParentElementOfDegree(polydeg2) ;
  coord   = mesh->GetCoordinates() ;

  // set polynomial degree to 'polydeg1' and 'polydeg2' in all elements
  coord->SetParentElements(parent1,1) ;
  coord->SetParentElements(parent2,2) ;

  // refresh the values
  coord->SetToCoordinates() ;

  //    2. Velocity field (GLL identical to coordinate field)
  //    -----------------------------------------------------

  velocity = new MortaredField(2,coord) ;
  
  //    3. Pressure field (GL)
  //    ----------------------

  shift          = 2 ;
  parents        = new Vector<ParentElement*>(2) ;
  polydegP       = max(polydeg1-shift,0) ;
  parents->At(1) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
  polydegP       = max(polydeg2-shift,0) ;
  parents->At(2) = ParentElementGL::GetParentElementOfDegree(polydegP) ;

  pressure = new FlatField(1,mesh,parents) ;

  delete parents ;
  

  // ====================
  // *** Forcing term ***
  // ====================

  FlatField *forcing ;
  Function  *forcingFunction ;

  forcing         = coord->DuplicateEmpty() ;
  forcingFunction = new ForcingNavierStokes2D() ;


  // ===========================
  // *** Boundary conditions ***
  // ===========================

  Mesh1D             *boundary ;
  FlatField          *prescribedV ;
  Function           *exactFunctionV, *exactFunctionP ;
  DirichletCondition *dcVx, *dcVy ;

  boundary = (Mesh1D*) mesh->GetSkin() ;
//  printf ("\nBoundary: ") ; boundary->Print() ;

  prescribedV = new FlatField(1,coord,boundary) ;
  dcVx        = new DirichletCondition(prescribedV,1) ;
  dcVy        = new DirichletCondition(prescribedV,2) ;
 
  exactFunctionV = new ExactStokes2D(1) ;          // 1 means "velocity"
  dcVx->SetFunction(exactFunctionV) ;
  dcVy->SetFunction(exactFunctionV) ;
  
  exactFunctionP = new ExactStokes2D(2) ;          // 2 means "pressure"


  // =========================
  // *** Integration rules ***
  // =========================

  FlatField *ruleMomentum, *ruleMass ;

  ruleMomentum = velocity->GetMain()->DuplicateEmpty(1) ;
  ruleMass     = pressure->DuplicateEmpty() ;


  // ======================================================
  // *** Linear solvers (overrides the default solvers) ***
  // ======================================================

  Solver *solverH, *solverP, *solverI ;

  solverH = new ConjugateGradient(500,1e-12) ;     // Helmholtz
  solverP = new ConjugateGradient(500,1e-12) ;     // pressure (external loop)
  solverI = new ConjugateGradient(500,1e-14) ;     // J=H-1 (case Uzawa only) 

  solverH->SetName("Helmholtz") ;
  solverP->SetName("pressure") ;
  solverI->SetName("internal") ;

  solverH->SetVerbosity(1) ;
  solverP->SetVerbosity(1) ;
  solverI->SetVerbosity(0) ;


  // =======================
  // *** Preconditioners ***
  // =======================

  Preconditioner *precondH, *precondP, *precondI ;

  precondH = NULL ;
//  precondH = new DiagonalPreconditioner() ;

  precondP = NULL ;
//  precondP = new DiagonalPreconditioner() ;
//  precondP = new MassPreconditioner() ;

  precondI = NULL ;
  precondI = new DiagonalPreconditioner() ;


  // ================================
  // *** Time-integration schemes ***
  // ================================

  TimeIntegrationScheme *diffusionScheme, *convectionScheme ;
 
# if defined BDF1
    diffusionScheme  = new BDF(1) ;
    convectionScheme = new Extrapolation(1) ;

# elif defined BDF2
    diffusionScheme  = new BDF(2) ;
    convectionScheme = new Extrapolation(2) ;

# elif defined BDF3
    diffusionScheme  = new BDF(3) ;
    convectionScheme = new Extrapolation(3) ;

# elif defined BDF4
    diffusionScheme  = new BDF(4) ;
    convectionScheme = new Extrapolation(4) ;

# endif

  diffusionScheme->SetTimeIncrement(dt) ;


  // =============================
  // *** Navier-Stokes problem ***
  // =============================

  NavierStokes *navierStokes ;
  SteadyStokes *steadyStokes ;
  FlatField    *startV0, *startV1, *startV2, *startV3 ;

  //    1. Main settings
  //    ----------------

  navierStokes = new NavierStokes(mesh,DECOMPOSITION_TYPE) ;
  navierStokes->SetTimeIntegrationScheme(diffusionScheme) ;
  navierStokes->SetConvectionScheme(convectionScheme) ;
  navierStokes->SetVelocity(velocity) ;
  navierStokes->SetPressure(pressure) ;
  navierStokes->SetReynolds(1.) ;
  navierStokes->SetForcing(forcing) ;
  navierStokes->SetIntegrationRuleMomentum(ruleMomentum) ;
  navierStokes->SetIntegrationRuleMass(ruleMass) ;
  navierStokes->SetForcingFunction(forcingFunction) ;
  navierStokes->AddDirichletConditionV(dcVx) ;
  navierStokes->AddDirichletConditionV(dcVy) ;

  //    2. Linear solvers and preconditioners
  //    -------------------------------------

  steadyStokes = navierStokes->GetSteadyStokes() ;

  // assign the linear solvers
  steadyStokes->GetHelmholtz()->SetSolver(solverH) ;
  steadyStokes->SetSolver(solverP) ;
  if (steadyStokes->GetInternalHelmholtz() != NULL)               // case Uzawa
    steadyStokes->GetInternalHelmholtz()->SetSolver(solverI) ;

  // assign the preconditioners
  steadyStokes->GetHelmholtz()->SetPreconditioner(precondH) ;
  steadyStokes->SetPreconditioner(precondP) ;
  if (steadyStokes->GetInternalHelmholtz() != NULL)              // case Uzawa
    steadyStokes->GetInternalHelmholtz()->SetPreconditioner(precondI) ;

  //    3. Initial values (and starting values, if multistep method)
  //    ------------------------------------------------------------

  startV0 = velocity->GetMain()->DuplicateEmpty() ;         // v(0)
  startV1 = velocity->GetMain()->DuplicateEmpty() ;         // v(-1)
  startV2 = velocity->GetMain()->DuplicateEmpty() ;         // v(-2)
  startV3 = velocity->GetMain()->DuplicateEmpty() ;         // v(-3)
  startV0->SetTo(exactFunctionV,0.) ;
  startV1->SetTo(exactFunctionV,-dt) ;
  startV2->SetTo(exactFunctionV,-dt-dt) ;
  startV3->SetTo(exactFunctionV,-dt-dt-dt) ;

# if defined BDF1
    navierStokes->SetValues(0,startV0) ;                    // v(0)

# elif defined BDF2
    navierStokes->SetValues(0,startV0) ;                    // v(0)
    navierStokes->SetValues(1,startV1) ;                    // v(-1)
    diffusionScheme->UseStarter(false) ;

# elif defined BDF3
    navierStokes->SetValues(0,startV0) ;                    // v(0)
    navierStokes->SetValues(1,startV1) ;                    // v(-1)
    navierStokes->SetValues(2,startV2) ;                    // v(-2)
    diffusionScheme->UseStarter(false) ;

# elif defined BDF4
    navierStokes->SetValues(0,startV0) ;                    // v(0)
    navierStokes->SetValues(1,startV1) ;                    // v(-1)
    navierStokes->SetValues(2,startV2) ;                    // v(-2)
    navierStokes->SetValues(3,startV3) ;                    // v(-3)
    diffusionScheme->UseStarter(false) ;

# else
    Error("main()","starting values not given") ;

# endif


  // ========================
  // *** Solution process ***
  // ========================

  FlatField *exactV, *exactP, *errorV, *errorP ;
  FILE      *converg ;
  real      t, normV, normP ;
  int       i ;
    
  exactV = velocity->GetMain()->DuplicateEmpty() ;
  exactP = pressure->DuplicateEmpty() ;
  errorV = exactV->DuplicateEmpty() ;
  errorP = exactP->DuplicateEmpty() ;

  for (i=1 ; i<=nsteps ; i++) {
    printf("\n\nStart step %d (time = %.4f)\n\n",
            i,diffusionScheme->GetTime()+dt) ;

    // approximate solution (v and p)

    navierStokes->Solve(1) ;                      // solve 1 step at a time

    if (i == nsteps) {

      // exact solution (v and p)

      t = diffusionScheme->GetTime() ;
      exactV->SetTo(exactFunctionV,t) ;
      exactP->SetTo(exactFunctionP,t) ;
      exactP->SetToZeroMean(1) ;

      // error on v and p

      errorV->CopyFrom(exactV) ;
      errorP->CopyFrom(exactP) ;
      errorV->Subtract(velocity->GetMain()) ;
      errorP->Subtract(pressure) ;
  
      // norm of error on v and p, and printings

      normV = errorV->GetH1Norm() / exactV->GetH1Norm() ;
      printf("\nH1norm(v-v_exact) / H1norm(v_exact) = %.4e",normV) ;

      normP = errorP->GetL2Norm() ;
      printf("\nL2norm(p-p_exact)                   = %.4e\n\n",normP) ;
    }

    if (i == nsteps) {

      // printings at last time step

      converg = Fopen("navier2DBP1PCBDF2.dat","a") ;
      fprintf(converg,"%.4e  %.4e %.4e  %d %d\n",
              dt,normV,normP,polydeg1,polydeg2) ;
      fclose(converg) ;
    }
  }

  // =============================================================
  // *** Optional postprocessing (Tecplot field visualization) ***
  // =============================================================

# ifdef POSTPROC

  Tecplot       *tecplot ;
  FlatField     *velocity2, *pressure2 ;
  ParentElement *parent ;

  // reinterpolate the fields on a fine grid to get nice plots (Tecplot uses
  // only straight lines)

  parent    = ParentElementGLL::GetParentElementOfDegree(20) ;
  velocity2 = new FlatField(2,mesh,parent) ;
  pressure2 = new FlatField(1,mesh,parent) ;

  velocity2->CopyInterpolateFrom(velocity->GetMain()) ;
  pressure2->CopyInterpolateFrom(pressure) ; 

  // reinterpolate similarly the coordinate field, otherwise Tecplot fails

  mesh->SetDefaultParentElement(parent) ;

  tecplot = new Tecplot("tecplot/navier2D_fields.dat",mesh->GetCoordinates()) ;
  tecplot->Put(velocity2,1,"vX") ;
  tecplot->Put(velocity2,2,"vY") ;
  tecplot->Put(pressure2,1,"Pressure") ;
  tecplot->GenerateOutput() ;

  delete tecplot ;
  unref(velocity2) ;
  unref(pressure2) ;

# endif

  // ================
  // *** Cleaning ***
  // ================

  Mpi_Finalize() ;

  return 0 ;
}

