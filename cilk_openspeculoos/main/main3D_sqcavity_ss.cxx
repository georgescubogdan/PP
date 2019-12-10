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

// main3D_sqcavity_ss.cxx
//
// Solves the lid-driven 3D cubic cavity for the steady Stokes equations.
//
//        0    = - grad p + nu Laplacian v
//     - div v = 0
//
// The data are:
//   . space domain: the unit square cavity ]0,1[ x ]0,1[ x ]0,1[;
//   . nu = 1;
//   . Dirichlet conditions:
//      - on the velocity: at all times, v=0 on all edges, except on the lid
//        (i.e., top face z=1), where it is prescribed in the X direction by a
//        regularized function:
//
//                      [ g1(x) g2(y) ]
//           v(x,y,1) = [      0      ]
//                      [      0      ]
//
//           g1(x) = (x*(1-x))**m1 * 2**(2 m1),
//           g2(y) = (y*(1-y))**m2 * 2**(2 m2),
//        where m1 = polydeg1/2, so g(x) is of degree polydeg1 (or polydeg1-1 
//        if polydeg1 is even), and similarly m2 = polydeg2/2.
//        Note: g(0) = g(1) = 0, g(0.5) = 1, for g=g1 and g=g2.
//      - on the pressure: none.


// Define the following option if you want the solution fields to be written
// into a file for Tecplot visualization:

//#define POSTPROC



#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //   2. nelem2:   number of elements in the Y direction.
  //   3. nelem3:   number of elements in the Z direction.
  //   4. polydeg1: polynomial degree of the GLL rule for the velocity field in
  //                the X direction, in all elements.
  //   5. polydeg2: polynomial degree of the GLL rule for the velocity field in
  //                the Y direction, in all elements.
  //   6. polydeg3: polynomial degree of the GLL rule for the velocity field in
  //                the Z direction, in all elements.

  int nelem1, nelem2, nelem3, polydeg1, polydeg2, polydeg3 ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem1   = 4 ; 
      nelem2   = 4 ; 
      nelem3   = 4 ; 
      polydeg1 = 5 ;
      polydeg2 = 5 ;
      polydeg3 = 5 ;
      break ;
    case 7 :
      nelem1   = atoi(argv[1]) ;
      nelem2   = atoi(argv[2]) ;
      nelem3   = atoi(argv[3]) ;
      polydeg1 = atoi(argv[4]) ;
      polydeg2 = atoi(argv[5]) ;
      polydeg3 = atoi(argv[6]) ;
      break ;
    default :
      printf("%s%s%s", "valid number of arguments:\n",
             "  0 or\n",
             "  6 (nelem1,nelem2,nelem3,polydeg1,polydeg2,polydeg3).\n\n") ;
      exit(0) ;
      break ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8 ;
  Vertex       *v1, *v2, *v3, *v4, *v5, *v6, *v7, *v8 ;
  Edge         *edge1, *edge2, *edge3, *edge4, *edge5, *edge6,
               *edge7, *edge8, *edge9, *edge10, *edge11, *edge12 ;
  Quad         *face1, *face2, *face3, *face4, *face5, *face6 ;
  Volume       *volume ;
  Distribution *distrib1, *distrib2, *distrib3 ;
  Mesh3D       *mesh ;

  //    1. Creation of the generation element
  //    -------------------------------------

  p1     = new Point(0., 0., 0.) ;
  p2     = new Point(1., 0., 0.) ;
  p3     = new Point(1., 1., 0.) ;
  p4     = new Point(0., 1., 0.) ;
  p5     = new Point(0., 0., 1.) ;
  p6     = new Point(1., 0., 1.) ;
  p7     = new Point(1., 1., 1.) ;
  p8     = new Point(0., 1., 1.) ;

  v1     = new Vertex(p1) ;
  v2     = new Vertex(p2) ;
  v3     = new Vertex(p3) ;
  v4     = new Vertex(p4) ;
  v5     = new Vertex(p5) ;
  v6     = new Vertex(p6) ;
  v7     = new Vertex(p7) ;
  v8     = new Vertex(p8) ;

  edge1  = new Line(v1,v2) ;
  edge2  = new Line(v2,v3) ;
  edge3  = new Line(v3,v4) ;
  edge4  = new Line(v4,v1) ;
  edge5  = new Line(v5,v6) ;
  edge6  = new Line(v6,v7) ;
  edge7  = new Line(v7,v8) ;
  edge8  = new Line(v8,v5) ;
  edge9  = new Line(v1,v5) ;
  edge10 = new Line(v2,v6) ;
  edge11 = new Line(v3,v7) ;
  edge12 = new Line(v4,v8) ;

  face1  = new Quad(edge4,edge8,edge12,edge9 ) ;
  face2  = new Quad(edge2,edge6,edge10,edge11) ;
  face3  = new Quad(edge1,edge5,edge9 ,edge10) ;
  face4  = new Quad(edge3,edge7,edge11,edge12) ;
  face5  = new Quad(edge1,edge3,edge2 ,edge4 ) ;
  face6  = new Quad(edge5,edge7,edge8 ,edge6 ) ;

  volume = new Hexa(face1,face2,face3,face4,face5,face6) ;

  //    2. Mesh generation
  //    ------------------

  distrib1 = new UniformDistribution(nelem1) ;
  distrib2 = new UniformDistribution(nelem2) ;
  distrib3 = new UniformDistribution(nelem3) ;

  volume->SetDistribution(1,distrib1) ;
  volume->SetDistribution(2,distrib2) ;
  volume->SetDistribution(3,distrib3) ;

  mesh = volume->GenerateMesh() ;
  mesh->DispatchElementsByBlocks(nelem1,nelem2,nelem3) ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement          *parent1, *parent2, *parent3 ;
  FlatField              *coord, *pressure ;
  MortaredField          *velocity ;
  Vector<ParentElement*> *parents ;
  int                    polydegP, shift ;

  //    1. Parent element definitions and coordinate field
  //    --------------------------------------------------

  parent1 = ParentElementGLL::GetParentElementOfDegree(polydeg1) ;
  parent2 = ParentElementGLL::GetParentElementOfDegree(polydeg2) ;
  parent3 = ParentElementGLL::GetParentElementOfDegree(polydeg3) ;
  coord   = mesh->GetCoordinates() ;

  // set polynomial degree to 'polydeg1' and 'polydeg2' in all elements
  coord->SetParentElements(parent1,1) ;
  coord->SetParentElements(parent2,2) ;
  coord->SetParentElements(parent3,3) ;

  // refresh the values
  coord->SetToCoordinates() ;

  // compute immediately invJacobian (this operation uses temporarily a lot of
  // memory space, and so should better be completed as early as possible)
  mesh->GetInvJacobian() ;

  //    2. Velocity field (GLL identical to coordinate field)
  //    -----------------------------------------------------

  velocity = new MortaredField(3,coord) ;
  
  //    3. Pressure field (GL)
  //    ----------------------

  shift          = 2 ;
  parents        = new Vector<ParentElement*>(3) ;
  polydegP       = max(polydeg1-shift,0) ;
  parents->At(1) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
  polydegP       = max(polydeg2-shift,0) ;
  parents->At(2) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
  polydegP       = max(polydeg3-shift,0) ;
  parents->At(3) = ParentElementGL::GetParentElementOfDegree(polydegP) ;

  pressure = new FlatField(1,mesh,parents) ;

  delete parents ;
  

  // ===========================
  // *** Boundary conditions ***
  // ===========================

  Mesh               *boundary, *lid, *walls ;
  FlatField          *prescribedVx, *prescribedVy, *prescribedVz,
                     *prescribedVxLid ; 
  Function           *functionLid ;
  DirichletCondition *dcVx, *dcVy, *dcVz, *dcVxLid ; 

  //    1. 1D boundary meshes
  //    ---------------------

  boundary = mesh->GetSkin() ;

  lid    = boundary->Extract(3,1.) ;                 // the lid (plane Z=1.)
  walls  = boundary->Except(lid) ;                   // the walls and the floor


  //    2. Boundary conditions on the velocity
  //    --------------------------------------

  // vX=0 on the walls, vX=regularized on the lid
  // vY=0 on the whole boundary (walls + lid)
  // vZ=0 on the whole boundary (walls + lid)

  prescribedVx    = new FlatField(1,coord,walls) ;
  prescribedVxLid = new FlatField(1,coord,lid) ;
  prescribedVy    = new FlatField(1,coord,boundary) ;
  prescribedVz    = new FlatField(1,coord,boundary) ;

  functionLid     = new RegularizedProfile(2,polydeg1/2,polydeg2/2) ;
  prescribedVxLid->SetTo(functionLid) ;

/*
  printf ("\nPrescribed on v_x (walls): ") ;    prescribedVx->Print() ;
  printf ("\nPrescribed on v_x (lid): ") ;      prescribedVxLid->Print() ;
  printf ("\nPrescribed on v_y (boundary): ") ; prescribedVy->Print() ;
  printf ("\nPrescribed on v_z (boundary): ") ; prescribedVz->Print() ;
*/

  dcVx    = new DirichletCondition(prescribedVx,1) ;
  dcVxLid = new DirichletCondition(prescribedVxLid,1) ;
  dcVy    = new DirichletCondition(prescribedVy,2) ;
  dcVz    = new DirichletCondition(prescribedVz,3) ;


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

  solverH = new ConjugateGradient(1000,1e-14) ;     // Helmholtz
  solverP = new ConjugateGradient(1000,1e-12) ;     // pressure (external loop)
  solverI = new ConjugateGradient(1000,1e-14) ;     // J = H-1

  solverH->SetName("Helmholtz") ;
  solverP->SetName("pressure") ;
  solverI->SetName("internal") ;

  solverH->SetVerbosity(1) ;
  solverP->SetVerbosity(3) ;
  solverI->SetVerbosity(1) ;


  // =======================
  // *** Preconditioners ***
  // =======================

  Preconditioner *precondH, *precondP, *precondI ;

  precondH = new DiagonalPreconditioner() ;
  precondH = NULL ;

  precondP = new DiagonalPreconditioner() ;
  precondP = new MassPreconditioner() ;
  precondP = new CouzyPreconditioner() ;
  precondP = NULL ;

  precondI = NULL ; 
  precondI = new DiagonalPreconditioner() ;


  // =============================
  // *** Steady Stokes problem ***
  // =============================

  SteadyStokes *problem ;
  real         nu ;

  //    1. Main settings
  //    ----------------

  nu = 1 ;

  problem = new SteadyStokes(mesh,UZAWA) ;
  problem->SetVelocity(velocity) ;
  problem->SetPressure(pressure) ;
  problem->SetNu(nu) ;
  problem->AddDirichletConditionV(dcVx) ;
  problem->AddDirichletConditionV(dcVxLid) ;
  problem->AddDirichletConditionV(dcVy) ;
  problem->AddDirichletConditionV(dcVz) ;
  problem->SetIntegrationRuleMomentum(ruleMomentum) ;
  problem->SetIntegrationRuleMass(ruleMass) ;

  //    2. Preconditioners of the 3 linear solvers
  //    ------------------------------------------

  // assign the linear solvers
  problem->GetHelmholtz()->SetSolver(solverH) ;
  problem->SetSolver(solverP) ;
  problem->GetInternalHelmholtz()->SetSolver(solverI) ;

  // assign the preconditioners
  problem->GetHelmholtz()->SetPreconditioner(precondH) ;
  problem->SetPreconditioner(precondP) ;
  problem->GetInternalHelmholtz()->SetPreconditioner(precondI) ;


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer("solution process") ;
  real  normV ;

  timer.Start() ;

  problem->Solve() ;

  timer.Stop() ;

  normV = velocity->GetMain()->GetH1Norm() ;
  Printf("\nH1norm(v_n) = %.4e\n",normV) ;


  // ============================================
  // *** Print timer and number of iterations ***
  // ============================================

  FILE *file ;

  if (Mpi_rank == 0) {
    timer.Print() ;
    file = Fopen("./timer_stokes3D_ss.dat","a+") ;
    fprintf(file,
      "\n3D sq. cavity ss %dx%dx%d elements, degree %dx%dx%d, %d procs:\n",
      nelem1,nelem2,nelem3,polydeg1,polydeg2,polydeg3,Mpi_nbProcessors) ;
    timer.PrintIn(file) ;
    fclose(file) ;
  }

  if (Mpi_rank == 0)
    solverI->Print() ;


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
  velocity2 = new FlatField(3,mesh,parent) ;
  pressure2 = new FlatField(1,mesh,parent) ;

  velocity2->CopyInterpolateFrom(velocity->GetMain()) ;
  pressure2->CopyInterpolateFrom(pressure) ; 

  // reinterpolate similarly the coordinate field, otherwise Tecplot fails

  mesh->SetDefaultParentElement(parent) ;

  tecplot = new Tecplot("tecplot/sqcavity2D_ss.dat",mesh->GetCoordinates()) ;
  tecplot->Put(velocity2,1,"vX") ;
  tecplot->Put(velocity2,2,"vY") ;
  tecplot->Put(velocity2,3,"vZ") ;
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
  unref(p5) ;
  unref(p6) ;
  unref(p7) ;
  unref(p8) ;
  unref(v1) ;
  unref(v2) ;
  unref(v3) ;
  unref(v4) ;
  unref(v5) ;
  unref(v6) ;
  unref(v7) ;
  unref(v8) ;
  unref(edge1) ;
  unref(edge2) ;
  unref(edge3) ;
  unref(edge4) ;
  unref(edge5) ;
  unref(edge6) ;
  unref(edge7) ;
  unref(edge8) ;
  unref(edge9) ;
  unref(edge10) ;
  unref(edge11) ;
  unref(edge12) ;
  unref(face1) ;
  unref(face2) ;
  unref(face3) ;
  unref(face4) ;
  unref(face5) ;
  unref(face6) ;
  unref(volume) ;
  unref(mesh) ;
  unref(distrib1) ;
  unref(distrib2) ;
  unref(distrib3) ;
  unref(velocity) ;
  unref(boundary) ;
  unref(lid) ;
  unref(walls) ;
  unref(prescribedVx) ;
  unref(prescribedVy) ;
  unref(prescribedVz) ;
  unref(prescribedVxLid) ;
  unref(prescribedVx) ;
  unref(dcVx) ;
  unref(dcVy) ;
  unref(dcVz) ;
  unref(ruleMomentum) ;
  unref(ruleMass) ;
  unref(dcVxLid) ;
  unref(solverH) ;
  unref(solverP) ;
  unref(solverI) ;
  unref(precondH) ;
  unref(precondI) ;
  unref(precondP) ;
  unref(problem) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}

