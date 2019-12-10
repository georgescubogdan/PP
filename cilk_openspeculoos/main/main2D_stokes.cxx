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

// main2D_stokes.cxx
//
//
// Solves one of the following 2D problems:
//
// - unsteady Stokes:
//
//       dv/dt = - grad p + nu Laplacian v + f(x,y,t)
//     - div v = 0
//
//     Exact solution: u(x,y,t) =   - cx sy st
//                     v(x,y,t) =     sx cy st
//                     p(x,y,t) = -pi sx sy st
//     with the notations:
//           cx = cos (pi/2 x), sx = sin (pi/2 x),
//           cy = cos (pi/2 y), sy = sin (pi/2 y),
//           ct = cos (pi/2 t), st = sin (pi/2 t).
//     Note: the mean pressure is 0.
//
// - steady Stokes:
//
//         0   = - grad p + nu Laplacian v + g(x,y)
//     - div v = 0
//
//     The exact solution is the same as that of the unsteady problem 
//     hereabove with t=1 (and thus st=1.).
//
// The data are:
//   . space domain: the square cavity ]0,1[ x ]0,1[ modelled by 2x2 elements;
//   . nu = 1;
//   . f(x,y,t) = [ - pi*pi cx sy st - pi/2 cx sy ct ] ;
//                [                    pi/2 sx cy ct ]
//   . g(x,y) = f(x,y,pi/2);
//   . Dirichlet conditions (on the velocity only): at all times, u and v are
//     set to the exact values on the 4 edges;
//   . initial conditions for the unsteady case: at t=0, v=0 and p=0 on the 
//     whole cavity (= exact values). 


// Define the following option if you want not only the steady, but also the
// unsteady case, to be computed:

//#define STEADY
#define UNSTEADY


// If you chose to run the unsteady case (ie, you defined UNSTEADY), define
// (i.e., uncomment) ONE among the following options (ignored if you did not 
// define UNSTEADY):

//#define BACKWARD_EULER
//#define FORWARD_EULER
//#define CRANK_NICOLSON
//#define BDF1                     // identical to BACKWARD_EULER
#define BDF2
//#define BDF3
//#define BDF4
//#define AB1                      // identical to FORWARD_EULER
//#define AB2
//#define AB3
//#define RK2
//#define RK4


// If you defined UNSTEADY and you use a non self-starting time-integration
// schemes (e.g., BDF2, BDF3, AB2), define the following option if you want to
// give the exact values for the starting step(s). If you do not define it, the
// default "starter" (typically, BDF1 or Crank-Nicolson) will be used.
// This option is ignored for self-starting time-integration schemes.

#define START_EXACT


// If you defined UNSTEADY, uncomment one among the following lines:
// (note: the steady problem can be solved only by the Uzawa method)

#define DECOMPOSITION_TYPE UZAWA
//#define DECOMPOSITION_TYPE BP1
//#define DECOMPOSITION_TYPE BP3
//#define DECOMPOSITION_TYPE BP1_PC
//#define DECOMPOSITION_TYPE BP3_PC


// Define the following option if you want the solution fields to be written
// into a file for Tecplot visualization:

#define POSTPROC



#include "core/main.hxx"



int main (int argc, char* argv[])
{ // The arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //   2. nelem2:   number of elements in the Y direction.
  //   3. polydeg1: polynomial degree of the GLL rule for the velocity field in
  //                the X direction, in all elements.
  //   4. polydeg2: polynomial degree of the GLL rule for the velocity field in
  //                the Y direction, in all elements.
  //   5. nsteps:   if you define option UNSTEADY: the number of time steps;
  //                if you undefine option UNSTEADY: ignored.
  //   6. dt:       if you define option UNSTEADY: the time increment delta_t;
  //                if you undefine option UNSTEADY: ignored.

  real dt ;
  int  nelem1, nelem2, polydeg1, polydeg2, nsteps ;

  Mpi_Initialize(&argc,&argv) ;

Printf("number of arguments: %d\n\n",argc) ;

  switch (argc) {
    case 1 :
      nelem1   = 4 ;
      nelem2   = 4 ;
      polydeg1 = 7 ;
      polydeg2 = 7 ;
      nsteps   = 20 ;
      dt       = 0.05 ;
      break ;
    case 5 :
      nelem1   = atoi(argv[1]) ;
      nelem2   = atoi(argv[2]) ;
      polydeg1 = atoi(argv[3]) ;
      polydeg2 = atoi(argv[4]) ;
      nsteps   = 1 ;
      dt       = 1. ;
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
      Printf("%s%s%s%s", "valid number of arguments:\n",
             "  0 or\n",
             "  4 (nelem1,nelem2,polydeg1,polydeg2) or\n\n",
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

  distrib1 = new UniformDistribution(nelem1) ;
  distrib2 = new UniformDistribution(nelem2) ;
  face->SetDistribution(1,distrib1) ;
  face->SetDistribution(2,distrib2) ;

  mesh = face->GenerateMesh() ;
  mesh->DispatchElementsByBlocks(nelem1,nelem2) ;

  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement          *parent1, *parent2 ;
  FlatField              *coord, *pressureS, *pressureU ;
  MortaredField          *velocityS, *velocityU ;
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

  velocityS = new MortaredField(2,coord) ;        // for the steady case
  velocityU = velocityS->Duplicate() ;            // for the unsteady case
  
  //    3. Pressure field (GL)
  //    ----------------------

  shift          = 2 ;
  parents        = new Vector<ParentElement*>(2) ;
  polydegP       = max(polydeg1-shift,0) ;
  parents->At(1) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
  polydegP       = max(polydeg2-shift,0) ;
  parents->At(2) = ParentElementGL::GetParentElementOfDegree(polydegP) ;

  pressureS = new FlatField(1,mesh,parents) ;     // for the steady case
  pressureU = pressureS->Duplicate() ;            // for the unsteady case


  // ====================
  // *** Forcing term ***
  // ====================

  FlatField *forcing ;
  Function  *forcingFunction ;

  forcing         = coord->DuplicateEmpty() ;
  forcingFunction = new ForcingStokes2D() ;


  // ===========================
  // *** Boundary conditions ***
  // ===========================

  Mesh1D             *boundary ;
  FlatField          *prescribedV ;
  Function           *exactFunctionV, *exactFunctionP ;
  DirichletCondition *dcVx, *dcVy ;

  boundary = (Mesh1D*) mesh->GetSkin() ;
  boundary->InterpolateCoordinatesLike(mesh) ;

  prescribedV = boundary->GetCoordinates()->DuplicateEmpty(1) ;
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

  ruleMomentum = velocityS->GetMain()->DuplicateEmpty(1) ;
  ruleMass     = pressureS->DuplicateEmpty() ;


  // =====================================================
  // *** Linear solvers (override the default solvers) ***
  // =====================================================

  Solver *solverH, *solverP, *solverI ;

  solverH = new ConjugateGradient(300,1e-14) ;     // Helmholtz
  solverP = new ConjugateGradient(5000,1e-12) ;    // pressure (external loop)
  solverI = new ConjugateGradient(300,1e-14) ;     // J=H-1 (case Uzawa only) 

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

  precondH = new DiagonalPreconditioner() ;
  //  precondH = NULL ;

  //  precondP = new MassPreconditioner() ;
  //  precondP = NULL ;
  precondP = new CouzyPreconditioner() ;
  //  precondP = new DiagonalPreconditioner() ;

  //  precondI = new MassPreconditioner() ;
  precondI = new DiagonalPreconditioner() ;
  //  precondI = NULL ;


  // ================
  // *** Problems ***
  // ================

  SteadyStokes          *steady ;
  Stokes                *unsteady ;
  TimeIntegrationScheme *scheme ;
  FlatField             *startV0, *startV1, *startV2, *startV3 ;
  real                  lambda, nu ;

  steady   = NULL ;
  unsteady = NULL ;
  scheme   = NULL ;
  startV0  = NULL ;
  startV1  = NULL ;
  startV2  = NULL ;
  startV3  = NULL ;

  lambda   = 0. ;
  nu       = 1. ;

  //    1. Steady Stokes
  //    ----------------

# ifdef STEADY

  steady = new SteadyStokes(mesh,UZAWA) ;
  steady->SetVelocity(velocityS) ;
  steady->SetPressure(pressureS) ;
  steady->SetLambda(lambda) ;
  steady->SetNu(nu) ;
  steady->SetForcing(forcing) ;
  steady->SetForcingFunction(forcingFunction) ;
  steady->AddDirichletConditionV(dcVx) ;
  steady->AddDirichletConditionV(dcVy) ;
  steady->SetIntegrationRuleMomentum(ruleMomentum) ;
  steady->SetIntegrationRuleMass(ruleMass) ;
  steady->SetSolver(solverP) ;
  steady->GetHelmholtz()->SetSolver(solverH) ;
  steady->GetInternalHelmholtz()->SetSolver(solverI) ;
  steady->SetPreconditioner(precondP) ;
  steady->GetHelmholtz()->SetPreconditioner(precondH) ;
  steady->GetInternalHelmholtz()->SetPreconditioner(precondI) ;

# endif

  //    2. Unsteady Stokes
  //    ------------------

# ifdef UNSTEADY

  // creation
  unsteady = new Stokes(mesh,DECOMPOSITION_TYPE) ;

  // time-integration scheme
# ifdef BACKWARD_EULER
    scheme = new ThetaMethod(1.) ;
# elif defined FORWARD_EULER
    scheme = new ThetaMethod(0.) ;
# elif defined CRANK_NICOLSON
    scheme = new ThetaMethod(0.5) ;
# elif defined BDF1
    scheme = new BDF(1) ;
# elif defined BDF2
    scheme = new BDF(2) ;
# elif defined BDF3
    scheme = new BDF(3) ;
# elif defined BDF4
    scheme = new BDF(4) ;
# elif defined AB1
    scheme = new AdamsBashforth(1) ;
# elif defined AB2
    scheme = new AdamsBashforth(2) ;
# elif defined AB3
    scheme = new AdamsBashforth(3) ;
# elif defined RK2
    scheme = new RungeKutta(2) ;
# elif defined RK4
    scheme = new RungeKutta(4) ;
# endif
  scheme->SetTimeIncrement(dt) ;

  // main settings
  unsteady->SetTimeIntegrationScheme(scheme) ;
  unsteady->SetVelocity(velocityU) ;
  unsteady->SetPressure(pressureU) ;
  unsteady->SetNu(nu) ;
  unsteady->SetForcing(forcing) ;
  unsteady->SetForcingFunction(forcingFunction) ;
  unsteady->AddDirichletConditionV(dcVx) ;
  unsteady->AddDirichletConditionV(dcVy) ;
  unsteady->SetIntegrationRuleMomentum(ruleMomentum) ;
  unsteady->SetIntegrationRuleMass(ruleMass) ;

  // linear solvers
  steady = unsteady->GetSteadyStokes() ;
  steady->GetHelmholtz()->SetSolver(solverH) ;
  steady->SetSolver(solverP) ;
  if (steady->GetInternalHelmholtz() != NULL)                   // case Uzawa
    steady->GetInternalHelmholtz()->SetSolver(solverI) ;

  // preconditioners
  steady->GetHelmholtz()->SetPreconditioner(precondH) ;
  steady->SetPreconditioner(precondP) ;
  if (steady->GetInternalHelmholtz() != NULL)                   // case Uzawa
    steady->GetInternalHelmholtz()->SetPreconditioner(precondI) ;

  // starting values
  startV0 = velocityU->GetMain()->DuplicateEmpty() ;
  startV1 = velocityU->GetMain()->DuplicateEmpty() ;
  startV2 = velocityU->GetMain()->DuplicateEmpty() ;
  startV3 = velocityU->GetMain()->DuplicateEmpty() ;
  startV0->SetTo(exactFunctionV,0.) ;                     // v(0)
  startV1->SetTo(exactFunctionV,-dt) ;                    // v(-1)
  startV2->SetTo(exactFunctionV,-dt-dt) ;                 // v(-2)
  startV3->SetTo(exactFunctionV,-dt-dt-dt) ;              // v(-3)

  // set v(0)
  scheme->SetValues(0,startV0) ;                          // v(0)

  // set v(-1), v(-2), etc, if desired
# ifdef START_EXACT  // use exact values for starting
# if defined BDF2
    scheme->SetValues(1,startV1) ;                        // v(-1)
    scheme->UseStarter(false) ;
# elif defined BDF3
    scheme->SetValues(1,startV1) ;                        // v(-1)
    scheme->SetValues(2,startV2) ;                        // v(-2)
    scheme->UseStarter(false) ;
# elif defined BDF4
    scheme->SetValues(1,startV1) ;                        // v(-1)
    scheme->SetValues(2,startV2) ;                        // v(-2)
    scheme->SetValues(3,startV3) ;                        // v(-3)
    scheme->UseStarter(false) ;
# endif
# endif   // end START_EXACT
   
# endif   // end UNSTEADY


  // ======================================
  // *** Solution process (steady case) ***
  // ======================================

# ifdef STEADY

  FlatField *exactV, *exactP, *errorV, *errorP ;
  FILE      *file ;
  real      normV, normP ;

  //    1. Solution
  //    -----------

  Printf("\nStart solving steady case\n\n") ;

  steady->SetTime(ONE) ;                           // steady solution at t=1.
  steady->Solve() ;

  //    2. Postprocessing
  //    -----------------

  // Check momemtum and mass equations

  MortaredField *sum = velocityS->DuplicateEmpty() ;
  FlatField     *Hv  = velocityS->GetMain()->DuplicateEmpty() ;
  FlatField     *Lv  = velocityS->GetMain()->DuplicateEmpty() ;
  FlatField     *Gp  = velocityS->GetMain()->DuplicateEmpty() ;
  FlatField     *Bf  = velocityS->GetMain()->DuplicateEmpty() ;
  FlatField     *Dv  = pressureS->DuplicateEmpty() ;

  Hv->SetToWeakIdentity(velocityS->GetMain(),ruleMomentum) ;
  Hv->Multiply(lambda) ;
  Lv->SetToWeakLaplacian(velocityS->GetMain(),ruleMomentum) ;
  Hv->Add(Lv,-nu) ;                // H = B - nu L

  Gp->SetToWeakGradient(pressureS,ruleMass) ;

  Bf->SetToWeakIdentity(forcing,ruleMomentum) ;

  sum->GetMain()->Add(Hv) ;        // sum = H v + G p - B f = 0 ?
  sum->GetMain()->Add(Gp) ;
  sum->GetMain()->Subtract(Bf) ;
  sum->Assemble() ;
  steady->GetHelmholtz()->ApplyHomogeneousDirichletConditions(sum);
  real normSum = sum->GetInfiniteNorm() ;
  Printf("\nnorm of H v + G p - B f: %.5e\n",normSum) ;
  
  Dv->SetToWeakDivergence(velocityS->GetMain(),ruleMass) ;
  real normDv = Dv->GetInfiniteNorm() ;
  Printf("\nnorm of Div v:           %.5e\n",normDv) ;

  pressureS->SetToZeroMean(1) ;

  exactV = velocityS->GetMain()->DuplicateEmpty() ;
  exactP = pressureS->DuplicateEmpty() ;
  exactV->SetTo(exactFunctionV,ONE) ;              // steady solution at t=1.
  exactP->SetTo(exactFunctionP,ONE) ;
  exactP->SetToZeroMean(1) ;
  
  errorV = exactV->Duplicate() ;
  errorP = exactP->Duplicate() ;
  errorV->Subtract(velocityS->GetMain()) ;
  errorP->Subtract(pressureS) ;

  normV = errorV->GetH1Norm() / exactV->GetH1Norm() ;
  Printf("\nH1norm(v-v_exact) / H1norm(v_exact) = %.4e\n",normV) ;

  normP = errorP->GetL2Norm() / exactP->GetL2Norm() ;
  Printf("\nL2norm(p-p_exact) / L2norm(p_exact) = %.4e\n\n",normP) ;

  file = Fopen("./stokes.dat","a") ;
  fprintf(file,"%d %.4e %.4e\n",polydeg1,normV,normP) ;
  fclose(file) ;

# endif     // end STEADY
  

  // ========================================
  // *** Solution process (unsteady case) ***
  // ========================================

# ifdef UNSTEADY

  FlatField *exactV, *exactP, *errorV, *errorP ;
  Timer     timer1("solution process") ;
  FILE      *converg ;
  real      t, normV, normP ;
  int       i ;
    
  exactV = velocityU->GetMain()->DuplicateEmpty() ;
  exactP = pressureU->DuplicateEmpty() ;
  errorV = exactV->Duplicate() ;
  errorP = exactP->Duplicate() ;

  Printf("\nStart solving unsteady case\n\n") ;

  for (i=1 ; i<=nsteps ; i++) {
    Printf("\n\nStart step %d (time = %.4f)\n",i,scheme->GetTime()+dt) ;

    // approximate solution (v and p)

    timer1.Start() ;
    scheme->DoOneStep() ;                            // solving process
    timer1.Stop() ;

    t = scheme->GetTime() ;

    if (i == nsteps) {

      // exact solution (v and p)

      exactV->SetTo(exactFunctionV,t) ;
      exactP->SetTo(exactFunctionP,t) ;

      // error on v and p

      errorV->CopyFrom(exactV) ;
      errorP->CopyFrom(exactP) ;
      errorP->SetToZeroMean(1) ;
      errorV->Subtract(velocityU->GetMain()) ;
      errorP->Subtract(pressureU) ;

      // norm of error on v and p, and printings

      normV = errorV->GetH1Norm() / exactV->GetH1Norm() ;
      Printf("\nH1norm(v-v_exact) / H1norm(v_exact) = %.4e\n",normV) ;

      normP = errorP->GetL2Norm() / exactP->GetL2Norm() ;
      Printf("\nL2norm(p-p_exact) / L2norm(p_exact) = %.4e\n\n",normP) ;
    }

    if (i == nsteps) {

      // printings at last time step

      if (DECOMPOSITION_TYPE == UZAWA)
        steady->GetInternalHelmholtz()->GetSolver()->Print() ;

      converg = Fopen("./stokes2DBP1PCBDF3.dat","a") ;
      fprintf(converg,"%.4e  %.4e %.4e  %d %d\n",
                       dt,normV,normP,polydeg1,polydeg2) ;
      fclose(converg) ;
    }
  }

  timer1.Print() ;

# endif    // end UNSTEADY


  // ======================
  // *** Postprocessing ***
  // ======================

# ifdef POSTPROC

  Tecplot       *tecplot ;
  FlatField     *velocity2, *pressure2 ;
  ParentElement *parent ;

  // reinterpolate the fields on a fine grid to get nice plots (Tecplot uses
  // only straight lines)

  parent    = ParentElementGLL::GetParentElementOfDegree(10) ;
  velocity2 = new FlatField(2,mesh,parent) ;
  pressure2 = new FlatField(1,mesh,parent) ;

  velocity2->CopyInterpolateFrom(steady->GetVelocity()->GetMain()) ;
  pressure2->CopyInterpolateFrom(steady->GetPressure()) ; 

  // reinterpolate similarly the coordinate field, otherwise Tecplot fails

  mesh->SetDefaultParentElement(parent) ;

  tecplot = new Tecplot("stokes2D.dat",mesh->GetCoordinates()) ;
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

  unref(p1) ;
  unref(p2) ;
  unref(p3) ;
  unref(p4) ;
  unref(v1) ;
  unref(v2) ;
  unref(v3) ;
  unref(v4) ;
  unref(edge1) ;
  unref(edge2) ;
  unref(edge3) ;
  unref(edge4) ;
  unref(face) ;
  unref(distrib1) ;
  unref(distrib2) ;
  unref(mesh) ;
  unref(velocityS) ;
  unref(velocityU) ;
  unref(pressureS) ;
  unref(pressureU) ;
  unref(forcing) ;
  unref(forcingFunction) ;
  unref(boundary) ;
  unref(prescribedV) ;
  unref(exactFunctionV) ;
  unref(exactFunctionP) ;
  unref(dcVx) ;
  unref(dcVy) ;
  unref(ruleMomentum) ;
  unref(ruleMass) ;
  unref(solverH) ;
  unref(solverP) ;
  unref(solverI) ;
  unref(precondH) ;
  unref(precondP) ;
  unref(precondI) ;
  unref(scheme) ;
  unref(steady) ;
  unref(unsteady) ;
  unref(startV0) ;
  unref(startV1) ;
  unref(startV2) ;
  unref(startV3) ;
  unref(exactV) ;
  unref(exactP) ;
  unref(errorV) ;
  unref(errorP) ;
  delete parents ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}
