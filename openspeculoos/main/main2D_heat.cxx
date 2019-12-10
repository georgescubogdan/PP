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

// main2D_heat.cxx
//
// Solves a 2D heat equation
//
//     dv/dt = Laplacian(v) + f(x,y,t)
//
// with:
//                     2      2    2      2
//     Laplacian(v) = d v / dx  + d v / dy
//
//                       3      -y
//     f(x,y,t) = pi/2 (x +3) (e  +1) cos(pi/2 t) 
//
//                     -y        3     -y
//             - (6x (e  +1) + (x +1) e  ) sin(pi/2 t)
//
//                                     3       -y
// The exact solution is: v(x,y,t) = (x + 3) (e  + 1) sin(pi/2 t).
//
// The data are:
//   . space domain: an irregular quadrilateral domain of which one edge is a
//     circular arch;
//   . Dirichlet conditions: v = v_exact on the 4 edges;
//   . initial conditions: at t=0, v=0 on the whole cavity.
//
// If the time-integration scheme is not self-starting (e.g., BDF3), the exact
// solution at starting steps (e.g., steps -1 and -2) are given.


// Define ONE among the following options, in order to choose the time-
// integration scheme to be used:

//#define BDF1
//#define BDF2
#define BDF3
//#define BDF4
//#define CRANK_NICOLSON


// Define (resp. undefine) the following option if you want (resp. do not want)
// the solution fields to be written into a file for Tecplot visualization:

//#define POSTPROC



#include "core/main.hxx"
#include <math.h>


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //   2. nelem2:   number of elements in the Y direction.
  //   3. polydeg1: polynomial degree of the GLL rule for v in the X direction,
  //                in all elements.
  //   4. polydeg2: polynomial degree of the GLL rule for v in the Y direction,
  //                in all elements.
  //   5. nsteps:   the number of time steps;
  //   6. dt:       time increment delta_t;

  real dt ;
  int  nelem1, nelem2, polydeg1, polydeg2, nsteps ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem1   = 24 ; 
      nelem2   = 16 ; 
      polydeg1 = 3 ;
      polydeg2 = 4 ;
      nsteps   = 50 ;
      dt       = 0.02 ;
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
  p3 = new Point(cos(PI/4.),cos(PI/4.)) ;
  p4 = new Point(0.,0.5) ;

  v1 = new Vertex(p1) ;
  v2 = new Vertex(p2) ;
  v3 = new Vertex(p3) ;
  v4 = new Vertex(p4) ;

  edge1 = new Line(v1,v2) ;
  edge2 = new Circle(v2,v3,p1) ;      // p1 = center of the circle
  edge3 = new Line(v3,v4) ;
  edge4 = new Line(v4,v1) ; 

  face  = new Quad(edge1,edge3,edge4,edge2) ;

  distrib1 = new UniformDistribution(nelem1) ;
  distrib2 = new UniformDistribution(nelem2) ;
  face->SetDistribution(1,distrib1) ;
  face->SetDistribution(2,distrib2) ;
  Quad::SetCurvedInteriorEdges(false) ;   // except on the circular edge, all

  mesh = face->GenerateMesh() ;

  mesh->DispatchElementsByBlocks(nelem1,nelem2) ;
  //  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent1, *parent2 ;
  FlatField     *coord ;
  MortaredField *v ;

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

  //    2. v (GLL identical to coordinate field)
  //    ---------------------------------------

  v = new MortaredField(1,coord) ;
  

  // ===========================
  // *** Boundary conditions ***
  // ===========================

  Mesh1D             *boundary ;
  FlatField          *prescribed ;
  Function           *exactFunction ;
  DirichletCondition *dc ;

  boundary = (Mesh1D*) mesh->GetSkin() ; 
  //  printf ("\nBoundary: ") ; boundary->Print() ;

  prescribed = new FlatField(1,coord,boundary) ;
  dc         = new DirichletCondition(prescribed,1) ;

  exactFunction = new ExactHeat2D() ;

  dc->SetFunction(exactFunction) ;


  // ====================
  // *** Forcing term ***
  // ====================

  FlatField *forcing ;
  Function  *forcingFunction ;

  forcing         = v->GetMain()->DuplicateEmpty() ;
  forcingFunction = new ForcingHeat2D() ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *rule ;

  rule = v->GetMain()->DuplicateEmpty(1) ;


  // =====================
  // *** Linear solver ***
  // =====================

  Solver *solver ;

  solver = new ConjugateGradient(1000,1e-14) ;

  solver->SetVerbosity(1) ;


  // ======================
  // *** Preconditioner ***
  // ======================

  Preconditioner *precond ;

  precond = NULL ;
  precond = new DiagonalPreconditioner() ;


  // ========================
  // *** Time-integration ***
  // ========================

  TimeIntegrationScheme *scheme ;

# if defined BDF1
    scheme = new BDF(1) ;

# elif defined BDF2
    scheme = new BDF(2) ;

# elif defined BDF3
    scheme = new BDF(3) ;

# elif defined BDF4
    scheme = new BDF(4) ;

# elif defined CRANK_NICOLSON
    scheme = new ThetaMethod(0.5) ;

# endif

  scheme->SetTimeIncrement(dt) ;


  // =============================
  // *** Heat-equation problem ***
  // =============================

  HeatEquation *problem ;
  real         alpha ;

  alpha = 1. ;

  problem = new HeatEquation(mesh) ;
  problem->SetField(v) ;
  problem->SetTimeIntegrationScheme(scheme) ;
  problem->SetThermalDiffusivity(alpha) ;
  problem->SetIntegrationRule(rule) ;
  problem->AddDirichletCondition(dc) ;
  problem->SetForcing(forcing) ;
  problem->SetForcingFunction(forcingFunction) ;
  problem->SetSolver(solver) ;
  problem->SetPreconditioner(precond) ;


  // =========================================
  // *** Starting values v(-1), v(-2), etc ***
  // =========================================

  FlatField             *startV0, *startV1, *startV2, *startV3 ;

  startV0 = v->GetMain()->DuplicateEmpty() ;
  startV1 = v->GetMain()->DuplicateEmpty() ;
  startV2 = v->GetMain()->DuplicateEmpty() ;
  startV3 = v->GetMain()->DuplicateEmpty() ;
  startV0->SetTo(exactFunction,0.) ;                     // v(0)
  startV1->SetTo(exactFunction,-dt) ;                    // v(-1)
  startV2->SetTo(exactFunction,-dt-dt) ;                 // v(-2)
  startV3->SetTo(exactFunction,-dt-dt-dt) ;              // v(-3)

  scheme->SetField(v) ;
  scheme->SetValues(0,startV0) ;                          // v(0)

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


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer("solution process") ;
  int   i ;

  timer.Start() ;

  for (i=1 ; i<=nsteps ; i++) {
    Printf("\n\nStart step %d (time = %.4f)\n\n",
            i,scheme->GetTime()+dt) ;

    problem->Solve(1) ;                           // solve 1 step at a time
  }

  timer.Stop() ;


  // ======================
  // *** Postprocessing ***
  // ======================

  FlatField *error ;
  FILE      *file ;
  real      norm ;

  //    1. Error on v
  //    -------------

  error = v->GetMain()->Duplicate() ;
//  printf("Solution: "); error->Print() ;

  error->SetTo(exactFunction,scheme->GetTime()) ;
//  printf("Exact solution: "); error->Print() ;

  error->Subtract(v->GetMain()) ;
//  printf("Error field: ") ; error->Print() ;

  norm = error->GetInfiniteNorm() ;
  Printf("\nInfinite norm of error field: %12.4E\n\n",norm) ;

  if (Mpi_rank == 0) {
    file = Fopen("./norm_heat_2D.dat","a+") ;
    fprintf(file,"%4d  %e\n",polydeg1,norm) ;
    fclose(file) ;
  }

  //    2. Print timer
  //    --------------

  if (Mpi_rank == 0) {
    timer.Print() ;
    file = Fopen("./timer_heat_2D.dat","a+") ;
    fprintf(file,"\n2D heat %dx%d elements, degree %dx%d, ",
                   nelem1,nelem2,polydeg1,polydeg2) ;
    fprintf(file,"%d steps, dt=%6.3f, %d procs\n",nsteps,dt,Mpi_nbProcessors);
    timer.PrintIn(file) ;
    fclose(file) ;
  }

  //   3. Optional postprocessing (Tecplot field visualization)
  //   --------------------------------------------------------

# ifdef POSTPROC

  Tecplot       *tecplot ;
  FlatField     *vBis ;
  ParentElement *parent ;

  // reinterpolate the fields on a fine grid to get nice plots (Tecplot uses
  // only straight lines)

  parent = ParentElementGLL::GetParentElementOfDegree(20) ;
  vBis   = new FlatField(2,mesh,parent) ;
  vBis->CopyInterpolateFrom(v->GetMain()) ;
Exit();

  // reinterpolate similarly the coordinate field, otherwise Tecplot fails

  mesh->SetDefaultParentElement(parent) ;

  tecplot = new Tecplot("heat2D_fields.dat",mesh->GetCoordinates()) ;
  tecplot->Put(vBis,1,"v") ;
  tecplot->GenerateOutput() ;
  delete tecplot ;

  unref(vBis) ;

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
  unref(v) ;
  unref(boundary) ;
  unref(prescribed) ;
  unref(exactFunction) ;
  unref(dc) ;
  unref(forcing) ;
  unref(forcingFunction) ;
  unref(rule) ;
  unref(solver) ;
  unref(precond) ;
  unref(scheme) ;
  unref(startV0) ;
  unref(startV1) ;
  unref(startV2) ;
  unref(startV3) ;
  unref(problem) ;
  unref(error) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}

