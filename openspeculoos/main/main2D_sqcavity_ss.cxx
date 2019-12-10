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

// main2D_sqcavity_ss.cxx
//
// Solves the lid-driven 2D square cavity for the steady Stokes equations.
//
//       0   = - grad p + nu Laplacian v
//    -div v = 0
//
// The data are:
//   . space domain: the unit square cavity ]0,1[ x ]0,1[;
//   . nu = 1;
//   . Dirichlet conditions:
//      - on the velocity: at all times, v=0 on all edges, except on the lid
//        (top edge), where it is prescribed by a regularized function
//           g(x) = (x*(1-x))**m * 2**2m,
//        where m = polydeg1/2, so g(x) is of degree polydeg1 (polydeg1-1 if
//        polydeg1 is even).
//        Note: g(0) = g(1) = 0, g(0.5) = 1.
//   . initial conditions: v=0 on the whole cavity.


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

  int nelem1, nelem2, polydeg1, polydeg2 ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem1   = 4 ; 
      nelem2   = 4 ; 
      polydeg1 = 7 ;
      polydeg2 = 7 ;
      break ;
    case 5 :
      nelem1   = atoi(argv[1]) ;
      nelem2   = atoi(argv[2]) ;
      polydeg1 = atoi(argv[3]) ;
      polydeg2 = atoi(argv[4]) ;
      break ;
    default :
      printf("%s%s%s", "valid number of arguments:\n",
             "  0 or\n",
             "  4 (nelem1,nelem2,polydeg1,polydeg2).\n\n") ;
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
  //    ----------------------------------------------------

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
  FlatField          *prescribedVx, *prescribedVxLid, *prescribedVy ;
  Function           *functionLid ;
  DirichletCondition *dcVx, *dcVxLid, *dcVy ; 
  Point              *midlid ;

  //    1. 1D boundary meshes
  //    ---------------------

  boundary = (Mesh1D*) mesh->GetSkin() ;

  midlid = new Point(0.5,1.) ;                       // the middle of the lid
  lid    = boundary->ExtractMeshBetween(p4,midlid,p3) ;
  walls  = boundary->ExtractMeshBetween(p4,p1,p3) ;  // the walls and the floor
/*
  printf ("\nBoundary: ") ; boundary->Print() ;
  printf ("\nLid: ") ; lid->Print() ;
  printf ("\nWalls: ") ; walls->Print() ;
*/

  //    2. Boundary conditions on the velocity
  //    --------------------------------------

  // vX=0 on the walls, vX=regularized on the lid
  // vY=0 on the whole boundary

  prescribedVx    = new FlatField(1,coord,walls) ;
  prescribedVxLid = new FlatField(1,coord,lid) ;
  prescribedVy    = new FlatField(1,coord,boundary) ;

  functionLid     = new RegularizedProfile(1,polydeg1/2) ;
  prescribedVxLid->SetTo(functionLid) ;

/*
  printf ("\nPrescribed on v_x (walls): ") ;    prescribedVx->Print() ;
  printf ("\nPrescribed on v_x (lid): ") ;      prescribedVxLid->Print() ;
  printf ("\nPrescribed on v_y (boundary): ") ; prescribedVy->Print() ;
*/

  dcVx    = new DirichletCondition(prescribedVx,1) ;
  dcVxLid = new DirichletCondition(prescribedVxLid,1) ;
  dcVy    = new DirichletCondition(prescribedVy,2) ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *ruleMomentum, *ruleMass ;

  ruleMomentum = velocity->GetMain()->DuplicateEmpty(1) ;
  ruleMass     = pressure->DuplicateEmpty() ;


  // ======================================================
  // *** Linear solvers (overrides the default solvers) ***
  // ======================================================

  Solver *solverH, *solverP, *solverI ;

  solverH = new ConjugateGradient(900,1e-12) ;     // Helmholtz
  solverP = new ConjugateGradient(900,1e-12) ;     // pressure (external loop)
  solverI = new ConjugateGradient(900,1e-14) ;     // J=H-1 (case Uzawa only) 

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
  precondP = new MassPreconditioner() ;

//  precondI = NULL ;
  precondI = new DiagonalPreconditioner() ;


  // =============================
  // *** Steady Stokes problem ***
  // =============================

  SteadyStokes *problem ;
  real         nu ;

  nu = 1. ;

  problem = new SteadyStokes(mesh,UZAWA) ;
  problem->SetVelocity(velocity) ;
  problem->SetPressure(pressure) ;
  problem->SetNu(nu) ;
  problem->AddDirichletConditionV(dcVx) ;
  problem->AddDirichletConditionV(dcVxLid) ;
  problem->AddDirichletConditionV(dcVy) ;
  problem->SetIntegrationRuleMomentum(ruleMomentum) ;
  problem->SetIntegrationRuleMass(ruleMass) ;
  problem->SetSolver(solverP) ;
  problem->GetHelmholtz()->SetSolver(solverH) ;
  problem->GetInternalHelmholtz()->SetSolver(solverI) ;
  problem->SetPreconditioner(precondP) ;
  problem->GetHelmholtz()->SetPreconditioner(precondH) ;
  problem->GetInternalHelmholtz()->SetPreconditioner(precondI) ;


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer("solution process") ;

  Printf("\n\nStart solving\n\n") ;
  timer.Start() ;

  problem->Solve() ;

  timer.Stop() ;

/*
  printf("\nSolution: velocity: ") ; steady->GetVelocity()->Print() ;
  printf("\nSolution: pressure: ") ; steady->GetPressure()->Print() ;
*/


  // ===================
  // *** Print timer ***
  // ===================

  FILE *file ;

  if (Mpi_rank == 0) {
    timer.Print() ;
    file = Fopen("./timer_ss_2D.dat","a+") ;
    fprintf(file,"\n2D sq. cavity %dx%d elements, %dx%d degree, %d procs:\n",
                  nelem1,nelem2,polydeg1,polydeg2,Mpi_nbProcessors) ;
    timer.PrintIn(file) ;
    fclose(file) ;
  }


  // ============================================================
  // *** Optional postprocessing: Tecplot field visualization ***
  // ============================================================

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

  tecplot = new Tecplot("tecplot/sqcavity2D_ss.dat",mesh->GetCoordinates()) ;
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
  unref(prescribedVx) ;
  unref(prescribedVxLid) ;
  unref(prescribedVy) ;
  unref(functionLid) ;
  unref(dcVx) ;
  unref(dcVxLid) ;
  unref(dcVy) ;
  unref(midlid) ;
  unref(ruleMomentum) ;
  unref(ruleMass) ;
  unref(solverH) ;
  unref(solverP) ;
  unref(solverI) ;
  unref(precondH) ;
  unref(precondP) ;
  unref(precondI) ;
  unref(problem) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}

