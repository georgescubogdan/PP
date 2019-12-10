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

// main3D_poisson.cxx
//
//
// Solves the 3D Poisson equation
//
//     u,   + u,   + u,   = f(x,y,z) ,
//       xx     yy     zz
//
//                      2
// with f(x,y,z) = -3 pi  cos (pi x) cos (pi y) cos(pi z)
//
// and Dirichlet conditions u(x,y,z) = u_exact(x,y,z) on the borders.
//
// The domain is a cube ]0,1[ x ]0,1[ x ]0,1[.
//
// The exact solution is: u_exact(x,y,z) = cos (pi x) cos (pi y) cos (pi z).


#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //   2. nelem2:   number of elements in the Y direction.
  //   3. nelem3:   number of elements in the Z direction.
  //   4. polydeg1: polynomial degree in the X direction.
  //   5. polydeg2: polynomial degree in the Y direction.
  //   6. polydeg3: polynomial degree in the Z direction.

  int nelem1, nelem2, nelem3, polydeg1, polydeg2, polydeg3 ;
  
  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem1   = 12 ;
      nelem2   = 8 ;
      nelem3   = 8 ;
      polydeg1 = 8 ;
      polydeg2 = 10 ;
      polydeg3 = 11 ;
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
      printf("0 argument or 6 arguments are required\n") ;
      exit(0) ;
      break ;
  }


  // =====================
  // *** Mesh creation ***
  // =====================

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

  Quad::SetCurvedInteriorEdges(false) ;     //   no transfinite edges

  mesh = volume->GenerateMesh() ;
  mesh->DispatchElementsByBlocks(nelem1,nelem2,nelem3) ;
//  mesh->PrintBarycenters() ;
//  mesh->PrintDispatching() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent1, *parent2, *parent3 ;
  FlatField     *coord ;
  MortaredField *vel ;

  //    1. Parent element definitions and coordinate field
  //    --------------------------------------------------

  parent1 = ParentElementGLL::GetParentElementOfDegree(polydeg1) ;
  parent2 = ParentElementGLL::GetParentElementOfDegree(polydeg2) ;
  parent3 = ParentElementGLL::GetParentElementOfDegree(polydeg3) ;
  coord = mesh->GetCoordinates() ;
  coord->SetParentElements(parent1,1) ;
  coord->SetParentElements(parent2,2) ;
  coord->SetParentElements(parent3,3) ;

  // refresh the values
  coord->SetToCoordinates() ;

  // compute immediately invJacobian (this operation uses temporarily a lot of
  // memory space, and so should better be completed as early as possible)
  mesh->GetInvJacobian() ;

  //    2. Field declarations and initializations
  //    -----------------------------------------

  vel = new MortaredField(1,coord) ;
  

  // ============================================
  // *** Forcing term and boundary conditions ***
  // ============================================

  Mesh               *boundary ;
  FlatField          *f, *prescribed ;
  Function           *functionForcing, *functionExact ;
  DirichletCondition *dc ; 

  //    1. Forcing term
  //    ---------------

  f               = vel->GetMain()->DuplicateEmpty(1) ;
  functionForcing = new Cos(-3.*PI*PI,PI,PI,PI) ;

  //    2. Boundary conditions (prescribed to the exact solution)
  //    ------------------------------------------------------------

  boundary = mesh->GetSkin() ;
//  printf ("\nBoundary: ") ; boundary->Print() ;

  prescribed    = boundary->GetCoordinates()->DuplicateEmpty(1) ;
  functionExact = new Cos(1.,PI,PI,PI) ;
  prescribed->SetTo(functionExact) ;

  dc = new DirichletCondition(prescribed,1) ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *rule ;

  rule = vel->GetMain()->DuplicateEmpty(1) ;


  // ==============
  // *** Solver ***
  // ==============

  Solver         *solver ;
  Preconditioner *precond ;

  solver = new ConjugateGradient(1000,1e-14) ;
  solver = new ConjugateGradient(1000,1e-12) ;
  solver->SetVerbosity(0) ;

  precond = new DiagonalPreconditioner() ;
  precond = NULL ;


  // ===============
  // *** Problem ***
  // ===============

  Helmholtz *problem ;

  problem = new Helmholtz(mesh) ;
  problem->SetAlpha(-1.) ;
  problem->SetUnknown(vel) ;
  problem->SetForcing(f) ;
  problem->SetForcingFunction(functionForcing) ;
  problem->AddDirichletCondition(dc) ;
  problem->SetIntegrationRule(rule) ;
  problem->SetSolver(solver) ;
  problem->SetPreconditioner(precond) ;


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer1("solution process") ;

  Printf("\n\nStart solving\n") ;
  timer1.Start() ;

  problem->Solve() ;

  timer1.Stop() ;

  FlatField::PrintNbFields() ;


  // ======================
  // *** Postprocessing ***
  // ======================

  FlatField *error ;
  FILE      *file ;
  real      H1Norm ;

  //    1. Timer
  //    --------

  if (Mpi_rank == 0) {
    timer1.Print() ;
    file = Fopen("./timer_poisson_3D.dat","a+") ;
    fprintf(file,"%dx%dx%d elements, degree %dx%dx%d, nb procs %d\n",
            nelem1,nelem2,nelem3,polydeg1,polydeg2,polydeg3,Mpi_nbProcessors) ;
    timer1.PrintIn(file) ;
    fclose(file) ;
  }

  //    2. Exact solution and error field
  //    ---------------------------------

  error = vel->GetMain()->GetWork(1) ;
  error->SetTo(functionExact) ;

  error->Subtract(vel->GetMain()) ;

  //    3. Error norm
  //    -------------

  H1Norm = error->GetH1Norm() ;
  if (Mpi_rank == 0) {
    printf("\nH1 norm of error field:        %12.4E\n\n",H1Norm) ;
    file = Fopen("./poisson3D.dat","a+") ;
    fprintf(file,"%d %e\n",polydeg1,H1Norm) ;
    fclose(file) ;
  }


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
  unref(vel) ;
  unref(f) ;
  unref(functionForcing) ;
  unref(boundary) ;
  unref(prescribed) ;
  unref(functionExact) ;
  unref(dc) ;
  unref(rule) ;
  unref(solver) ;
  unref(precond) ;
  error->GetOwner()->Retrieve(error) ;
  unref(problem) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}
