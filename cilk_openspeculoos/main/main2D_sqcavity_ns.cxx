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

// main2D_sqcavity_ns.cxx
//
// Solves the lid-driven 2D square cavity for the Navier-Stokes equations.
//
//     dv/dt + v.grad v = - grad p + 1/Re Laplacian v
//     - div v          = 0
//
// The data are:
//   . space domain: the unit square cavity ]0,1[ x ]0,1[;
//   . Re = 1000;
//   . Dirichlet conditions:
//      - on the velocity: at all times, v=0 on all edges, except on the lid
//        (top edge), where it is prescribed by a regularized function 
//        g(x) = 16 (x*(1-x))**2m.
//        Note: g(0) = g(1) = 0, g(0.5) = 1.
//      - on the pressure: none.
//   . initial conditions: at t=0, v=0 on the whole cavity.
//
// Karniadakis's method BDFn/EXn (n=1, 2, 3 or 4) is used. See Speculoos
// programmer's manual or Couzy's thesis, p 40.


// In order to choose the time-integration scheme to be used for the diffusive
// part of the problem, define (i.e., uncomment) ONE among the following
// options:

//#define BDF1
#define BDF2
//#define BDF3
//#define BDF4


// Define (i.e., uncomment) ONE among the following lines (the most accurate
// is the Uzawa one):

//#define DECOMPOSITION_TYPE UZAWA
//#define DECOMPOSITION_TYPE BP1
//#define DECOMPOSITION_TYPE BP3
#define DECOMPOSITION_TYPE BP1_PC


// Define the following option if you want the solution field (v,p) at the last
// time step to be written into a file for Tecplot visualization:

#define TECPLOT


// Define the following option if you want the stream function psi to be
// computed at the last time step. This requires solving a Poisson problem
//   Laplacian psi = -omega, where omega = rotational v,
// with homogeneous Dirichlet conditions.

//# define STREAM_FUNCTION




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
      cout << "Using default arguments" << endl;
      nelem1   = 8 ; 
      nelem2   = 8 ; 
      polydeg1 = 8 ;
      polydeg2 = 8 ;
      dt       = 0.001;
      nsteps   = (int)(1./dt)*3; // Three macroscopic steps.
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
  real reynolds = 100.;

  if (nelem1%2 != 0)
    Error("main()","'nelem1' must be even") ;


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
//  printf("\nCoordinates: ") ; coord->Print() ;

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
  

  // ===========================
  // *** Boundary conditions ***
  // ===========================

  Mesh1D             *boundary, *lid, *walls ;
  FlatField          *prescribedVxWalls, *prescribedVxLid, *prescribedVy ;
  Function           *functionLid ;
  DirichletCondition *dcVxWalls, *dcVxLid, *dcVy ; 
  Point              *midlid ;

  //    1. 1D boundary meshes
  //    ---------------------

  boundary = (Mesh1D*) mesh->GetSkin() ;

  //lid = (Mesh1D*) boundary->Extract(2,1.);
  //walls = (Mesh1D*) boundary->Except(lid);

  midlid = new Point(0.5,1.) ;                       // the middle of the lid
  lid    = boundary->ExtractMeshBetween(p4,midlid,p3) ;
  walls  = boundary->ExtractMeshBetween(p4,p1,p3) ;  // the walls and the floor


  //    2. Boundary conditions on the velocity
  //    --------------------------------------

  // vX=0 on the walls, vX=regularized on the lid
  // vY=0 on the whole boundary

  prescribedVxWalls = walls   ->GetCoordinates()->DuplicateEmpty(1) ;
  prescribedVxLid   = lid     ->GetCoordinates()->DuplicateEmpty(1) ;
  prescribedVy      = boundary->GetCoordinates()->DuplicateEmpty(1) ;

  functionLid       = new RegularizedProfile(1,18,1);
  prescribedVxLid->SetTo(functionLid) ;

  dcVxWalls = new DirichletCondition(prescribedVxWalls,1) ;
  dcVxLid   = new DirichletCondition(prescribedVxLid,1) ;
  dcVy      = new DirichletCondition(prescribedVy,2) ;


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
  precondH = new DiagonalPreconditioner() ;

  precondP = NULL ;
  //precondP = new MassPreconditioner() ;
  //precondP = new DiagonalPreconditioner() ;
  precondP = new CouzyPreconditioner() ;

  precondI = NULL ;
  precondI = new DiagonalPreconditioner() ;


  // ========================
  // *** Time-integration ***
  // ========================

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
  navierStokes->SetReynolds(reynolds) ;
  navierStokes->SetIntegrationRuleMass(ruleMass) ;
  navierStokes->SetIntegrationRuleMomentum(ruleMomentum) ;
  navierStokes->AddDirichletConditionV(dcVxWalls) ;
  navierStokes->AddDirichletConditionV(dcVxLid) ;
  navierStokes->AddDirichletConditionV(dcVy) ;

  steadyStokes = navierStokes->GetSteadyStokes() ;

  //    2. Preconditioners of the 3 linear solvers
  //    ------------------------------------------

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


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer("solution process") ;
  real  normV ;
  int   i ;
 
  timer.Start() ;
 
  for (i=1 ; i<=nsteps ; i++) {
    Printf("\n\nStart step %d (time = %.4f)\n\n",
            i,diffusionScheme->GetTime()+dt) ;

    navierStokes->Solve(1) ;                      // solve 1 step at a time
  }

  timer.Stop() ;

  normV = velocity->GetMain()->GetH1Norm() ;
  Printf("\nH1norm(v_n) = %.4e\n",normV) ;

  // ===================
  // *** Print timer ***
  // ===================

  FILE *file ;

  if (Mpi_rank == 0) {
    timer.Print() ;
    file = Fopen("./timer_sqcavity_ns.dat","a+") ;
    fprintf(file,
      "\n2D sq. cavity ns %dx%d elements, degree %dx%d, Re %5.1f, %d procs:\n",
      nelem1,nelem2,polydeg1,polydeg2,reynolds,Mpi_nbProcessors) ;
    timer.PrintIn(file) ;
    fclose(file) ;
  }


  // ============================================
  // *** Optional postprocessing: stream function
  // ============================================

# ifdef STREAM_FUNCTION

  Helmholtz          *poisson ;
  Solver             *solverPsi ;
  MortaredField      *stream ;
  FlatField          *omega ;
  DirichletCondition *dc ;

  stream = velocity->DuplicateEmpty(1) ;

  omega = velocity->GetMain()->DuplicateEmpty(1) ;
  omega->SetToRotational(velocity->GetMain()) ;

  dc = new DirichletCondition(prescribedVy,1) ;

  solverPsi = new ConjugateGradient(2000,1.e-12) ;
  solverPsi->SetVerbosity(3) ;

  poisson = new Helmholtz(mesh) ;
  poisson->SetAlpha(1.) ;              // to solve 0 = Laplacian stream + omega
  poisson->SetUnknown(stream) ;
  poisson->SetIntegrationRule(ruleMomentum) ;
  poisson->SetForcing(omega) ;
  poisson->AddDirichletCondition(dc) ;
  poisson->SetSolver(solverPsi) ;

  Printf("\nStart solving for stream function\n") ;

  poisson->Solve() ;                   // yields 'stream'

# else

  MortaredField *stream  = velocity->DuplicateEmpty(1) ;

# endif


  // ============================================================
  // *** Optional postprocessing: Tecplot field visualization ***
  // ============================================================

# ifdef TECPLOT

  Tecplot       *tecplot ;
  FlatField     *velocity2, *pressure2;
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

  // dump Tecplot file

  tecplot = new Tecplot("tecplot/sqcavity2D_ns.dat",mesh->GetCoordinates()) ;
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
  unref(velocity) ;
  unref(pressure) ;
  unref(boundary) ;
  unref(lid) ;
  unref(walls) ;
  unref(prescribedVxWalls) ;
  unref(prescribedVxLid) ;
  unref(prescribedVy) ;
  unref(functionLid) ;
  unref(dcVxWalls) ;
  unref(dcVxLid) ;
  unref(dcVy) ;
  unref(ruleMomentum) ;
  unref(ruleMass) ;
  unref(midlid) ;
  unref(solverH) ;
  unref(solverP) ;
  unref(solverI) ;
  unref(precondH) ;
  unref(precondP) ;
  unref(precondI) ;
  unref(diffusionScheme) ;
  unref(convectionScheme) ;
  unref(navierStokes) ;
  unref(stream) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}

