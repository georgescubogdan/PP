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

// main2D_poisson.cxx
//
// Solves the following 2D Poisson equation:
//
//    0 = - Laplacian u + f(x,y)
//
// with:
//    . Laplacian u = u,   + u,
//                      xx     yy
//
//    . f(x,y) = -2 pi cos (pi x) cos (pi y)
//
// on the rectangular domain [-1.0, 1.2] x [-1.3, 1.0].
//
// The exact solution is: u_exact(x,y) = cos (pi x) cos (pi y).
//
// The boundary conditions are:
//
//   . Dirichlet conditions on half of the contour (i.e., on 2 edges), with:
//       u(x,y) = u_exact(x,y) on the borders.
//
//   . Neumann conditions on the other half of the contour, with:
//       grad u(x,y) = [ du_exact/dx ] = [-pi sin(pi x) cos(pi y) ]
//                     [ du_exact/dy ] = [-pi cos(pi x) sin(pi y) ]
//
// The domain is chosen different to the more usual [-1,1] x [-1,1] in order to
// make the effect of Neumann conditions less trivial. Note that if you modify 
// the domain, you will also converge to the exact solution.


#include "core/main.hxx"


int main (int argc, char** argv)
{ // The (optional) arguments are:
  //   1. nelem1:   number of elements in the X direction.
  //                default = 1.
  //   2. nelem2:   number of elements in the Y direction.
  //                default = 1.
  //   3. polydeg1: polynomial degree of the GLL rule in the X direction.
  //                default = 2.
  //   4. polydeg2: polynomial degree of the GLL rule in the Y direction.
  //                default = 2.

  int nelem1, nelem2, polydeg1, polydeg2 ;

  Mpi_Initialize(&argc,&argv) ;
  
  switch (argc) {
    case 1 :
      nelem1   = 1 ;
      nelem2   = 1 ;
      polydeg1 = 2 ;
      polydeg2 = 2 ;
      break ;
    case 5 :
      nelem1   = atoi(argv[1]) ;
      nelem2   = atoi(argv[2]) ;
      polydeg1 = atoi(argv[3]) ;
      polydeg2 = atoi(argv[4]) ;
      break ;
    default :
      Error("main","0 or 4 arguments (nx,ny,polydeg1,polydeg2) are required") ;
      break ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2, *p3, *p4 ;
  Vertex       *v1, *v2, *v3, *v4 ;
  Edge         *line1, *line2, *line3, *line4 ;
  Face         *face ;
  Distribution *distrib1, *distrib2 ;
  Mesh2D       *mesh ;

  p1 = new Point(-1.0, -1.2) ;
  p2 = new Point( 1.3, -1.2) ;
  p3 = new Point( 1.3,  1.0) ;
  p4 = new Point(-1.0,  1.0) ;

  v1 = new Vertex(p1) ;
  v2 = new Vertex(p2) ;
  v3 = new Vertex(p3) ;
  v4 = new Vertex(p4) ;

  line1 = new Line(v1,v2) ;
  line2 = new Line(v2,v3) ;
  line3 = new Line(v3,v4) ;
  line4 = new Line(v1,v4) ;

  face  = new Quad(line1,line3,line4,line2) ;

  distrib1 = new UniformDistribution(nelem1) ;
  distrib2 = new UniformDistribution(nelem2) ;
//distrib1 = new GeometricalDistribution(nelem1,0.40) ;
//distrib2 = new GeometricalDistribution(nelem2,0.45) ;
  face->SetDistribution(1,distrib1) ;
  face->SetDistribution(2,distrib2) ;

  mesh = face->GenerateMesh() ;
  mesh->DispatchElementsByBlocks(nelem1,nelem2) ;
  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent1, *parent2 ;
  FlatField     *coord ;
  MortaredField *vel ;

  //    1. Parent element definitions and coordinate field
  //    --------------------------------------------------

  parent1 = ParentElementGLL::GetParentElementOfDegree(polydeg1) ;
  parent2 = ParentElementGLL::GetParentElementOfDegree(polydeg2) ;
  coord = mesh->GetCoordinates() ;
  coord->SetParentElements(parent1,1) ;
  coord->SetParentElements(parent2,2) ;
  coord->SetToCoordinates() ;                        // refresh the values
//  printf("\nCoordinates :") ; coord->Print() ;

  //    2. Field declarations and initializations
  //    -----------------------------------------

  vel = new MortaredField(1,coord) ;
  

  // ============================================
  // *** Forcing term and boundary conditions ***
  // ============================================

  Mesh1D             *contour, *dirichletMesh, *neumannMesh ;
  FlatField          *f, *prescribed, *prescribedGradient, *ruleContour ;
  Function           *forcingFunction, *exactFunction, *exactGradientXFunction,
                     *exactGradientYFunction ;
  DirichletCondition *dc ; 
  NeumannCondition   *nc ;

  //    1. Forcing term
  //    ---------------

  f = vel->GetMain()->DuplicateEmpty() ;

  forcingFunction = new Cos(-2.*PI*PI,PI,PI) ;
  f->SetTo(forcingFunction) ;

  //    2. Dirichlet conditions (prescribed to the exact solution)
  //    ---------------------------------------------------------

  contour = (Mesh1D*) mesh->GetSkin() ;

  dirichletMesh = contour->ExtractMeshBetween(p1,p4,p3) ; 

  prescribed    = dirichletMesh->GetCoordinates()->DuplicateEmpty(1) ;

  exactFunction = new Cos(1.,PI,PI) ;
  prescribed->SetTo(exactFunction) ;

  dc = new DirichletCondition(prescribed,1) ;

  //    3. Neumann conditions (prescribed to the exact solution gradient)
  //    -----------------------------------------------------------------

  neumannMesh            = contour->ExtractMeshBetween(p1,p2,p3) ;

  prescribedGradient     = neumannMesh->GetCoordinates()->DuplicateEmpty(2) ;

  exactGradientXFunction = new SinCos(-PI,PI,PI) ;
  exactGradientYFunction = new CosSin(-PI,PI,PI) ;

  prescribedGradient->SetTo(exactGradientXFunction,1,1) ;
  prescribedGradient->SetTo(exactGradientYFunction,2,1) ;

  nc = new NeumannCondition(prescribedGradient,1,ONE,ONE_MINUS) ;

  ruleContour = neumannMesh->GetCoordinates()->DuplicateEmpty(1) ;
  nc->SetIntegrationRule(ruleContour) ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *rule ;

  rule = vel->GetMain()->DuplicateEmpty(1) ;


  // =================================
  // *** Solver and preconditioner ***
  // =================================

  Solver *solver ;
  Preconditioner *precond ;

  solver = new ConjugateGradient(1000,1.e-14) ;
  solver->SetVerbosity(1) ;

  precond = new MassPreconditioner() ;
  precond = new EuclidianNormPreconditioner() ;
  precond = NULL ;
  precond = new DiagonalPreconditioner() ;


  // ===============
  // *** Problem ***
  // ===============

  Helmholtz *problem ;

  problem = new Helmholtz(mesh) ;
  problem->SetAlpha(-1.) ;
  problem->SetUnknown(vel) ;
  problem->SetForcing(f) ;
  problem->SetIntegrationRule(rule) ;
  problem->AddDirichletCondition(dc) ;
  problem->AddNeumannCondition(nc) ;
  problem->SetSolver(solver) ;
  problem->SetPreconditioner(precond) ;


  // ========================
  // *** Solution process ***
  // ========================

  Timer timer1("solution process") ;

  timer1.Start() ;
  problem->Solve() ;
  timer1.Stop() ;


  // ======================
  // *** Postprocessing ***
  // ======================

  FlatField *error ;
  real      infNorm, h1Norm ;

  error = vel->GetMain()->DuplicateEmpty() ;

  error->SetTo(exactFunction) ;
//  printf("Exact solution: "); error->Print() ;

  error->Subtract(vel->GetMain()) ;
  //  printf("Error field: ") ; error->Print() ;

  infNorm = error->GetInfiniteNorm() ;
  h1Norm  = error->GetH1Norm() ;
  Printf("\nInfinite norm of error field: %12.4E",infNorm) ;
  Printf("\nH1 norm of error field:       %12.4E\n\n",h1Norm) ;

  if (Mpi_rank == 0)
    timer1.Print() ;


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
  unref(line1) ;
  unref(line2) ;
  unref(line3) ;
  unref(line4) ;
  unref(distrib1) ;
  unref(distrib2) ;
  unref(face) ;
  unref(mesh) ;
  unref(vel) ;
  unref(f) ;
  unref(forcingFunction) ;
  unref(contour) ;
  unref(dirichletMesh) ;
  unref(prescribed) ;
  unref(exactFunction) ;
  unref(dc) ;
  unref(neumannMesh) ;
  unref(prescribedGradient) ;
  unref(exactGradientXFunction) ;
  unref(exactGradientYFunction) ;
  unref(nc) ;
  unref(ruleContour) ;
  unref(rule) ;
  unref(problem) ;
  unref(precond) ;
  unref(solver) ;
  unref(error) ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}
