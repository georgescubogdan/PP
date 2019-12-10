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

// main1D_poisson.cxx
//
//
// Solves the following 1D Poisson equation:
//
//                                  2
//     u,   = f(x) , with f(x) = -pi cos(pi x)
//       xx
//
//     and boundary conditions u(-1) = u(1) = -1
//     on a mesh [-1,1].
//     The exact solution is: u(x) = cos(pi x).


// Define (resp. undefine) USE_CG if you want a conjugate-gradient solver
// (resp. a direct solver) for solving the linear system. The output is the
// same. The standard solver is the C.G. one; however, by using a direct 
// solver the user can figure out more easily how the solution phase works.

#define USE_CG


#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. number of elements 
  //   2. polynomial degree of GLL rule.

  int nelem, polydeg ;

  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem   = 1 ;                        // default value
      polydeg = 2 ;                        // default value
      break ;
    case 2 :
      nelem   = atoi(argv[1]) ;
      polydeg = 2 ;                        // default value
      break ;
    case 3 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      break ;
    default :
      Error("main","at most 2 arguments are accepted") ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2 ;
  Vertex       *v1, *v2 ;
  Edge         *line ;
  Distribution *distrib ;
  Mesh1D       *mesh ;

  p1 = new Point(-1.) ;
  p2 = new Point( 1.) ;

  v1 = new Vertex(p1) ;
  v2 = new Vertex(p2) ;

  line = new Line(v1,v2) ;

  distrib = new UniformDistribution(nelem) ;
  line->SetDistribution(1,distrib) ;

  mesh = line->GenerateMesh() ;
  mesh->DispatchElements() ;
  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent ;
  FlatField     *coord ;
  MortaredField *vel ;

  //    1. Parent element definitions
  //    -----------------------------

  parent = ParentElementGLL::GetParentElementOfDegree(polydeg) ;
  mesh->SetDefaultParentElement(parent) ;

  //    2. Field declarations and initializations
  //    -----------------------------------------

  coord = mesh->GetCoordinates() ;
  vel = new MortaredField(1,coord) ;


  // ============================================
  // *** Forcing term and boundary conditions ***
  // ============================================

  Mesh               *boundary ;
  FlatField          *f, *prescribed ;
  Function           *function1, *function2 ;
  DirichletCondition *dc ; 

  //    1. Forcing term
  //    ---------------

  f         = vel->GetMain()->DuplicateEmpty() ;
  function1 = new Cos(-PI*PI,PI) ;
  f->SetTo(function1) ;
  
  //    2. Boundary conditions
  //    ----------------------

  boundary = mesh->GetSkin() ;
//  printf ("\nBoundary: ") ; boundary->PrintBarycenters() ;

  prescribed = boundary->GetCoordinates()->DuplicateEmpty() ;
  function2  = new Cos(ONE,PI) ;
  prescribed->SetTo(function2) ;
//  printf ("\nPrescribed: "); prescribed->Print() ;

  dc = new DirichletCondition(prescribed,1) ;


  // ========================
  // *** Integration rule ***
  // ========================

  FlatField *rule ;

  rule = vel->GetMain()->DuplicateEmpty(1) ;


  // =================================
  // *** Solver and preconditioner ***
  // =================================

  Solver         *solver ;
  Preconditioner *precond ;

# ifdef USE_CG
    solver = new ConjugateGradient(100,1e-12) ;
    solver->SetVerbosity(3) ;
# else
    solver = new DirectSolver() ;
# endif

    //  precond = NULL ;
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
  problem->SetSolver(solver) ;
  problem->SetPreconditioner(precond) ;


  // ========================
  // *** Solution process ***
  // ========================

  problem->Solve() ;

//  printf("\nSolution: ") ; problem->GetUnknown()->Print() ;


  // ======================
  // *** Postprocessing ***
  // ======================

  FlatField *exact, *error ;
  FILE      *tecfile ;
  Tecplot   *postProcessor ;
  real      infiniteNorm, H1Norm ;

  //    1. Exact solution and error field
  //    ---------------------------------

  exact = vel->GetMain()->DuplicateEmpty(1) ;
  exact->SetTo(function2) ;

  error = exact->Duplicate() ;
  error->Subtract(vel->GetMain()) ;
//  printf("\nError field: ") ; error->Print() ;

  //    2. Error norm
  //    -------------

  infiniteNorm = error->GetInfiniteNorm() ;
  H1Norm       = error->GetH1Norm() ;

  printf("\nInfinite norm of error field: %12.4E",infiniteNorm) ;
  printf("\nH1 norm of error field:       %12.4E\n\n",H1Norm) ;

  if (Mpi_rank == 0) {
    tecfile = Fopen("poisson1D_norm.dat","a+") ;
    fprintf(tecfile,"%d %e\n",polydeg,H1Norm) ;
    fclose(tecfile) ;
  }

  //    3. Tecplot data
  //    ---------------

  if (Mpi_nbProcessors == 1) {
    postProcessor = new Tecplot("poisson1D_files.dat",coord) ;
    postProcessor->SetTitle("1D Poisson") ;  
    postProcessor->Put(vel->GetMain(),1,"velocity") ;  
    postProcessor->Put(exact,1,"exact_velocity") ;
    postProcessor->Put(error,1,"error") ;
    postProcessor->GenerateOutput() ;
  }
  else
    postProcessor = NULL ;


  // ================
  // *** Cleaning ***
  // ================

  unref(p1) ;
  unref(p2) ;
  unref(v1) ;
  unref(v2) ;
  unref(line) ;
  unref(distrib) ;
  unref(mesh) ;
  unref(vel) ;
  unref(f) ;
  unref(prescribed) ;
  unref(boundary) ;
  unref(dc) ;
  unref(rule) ;
  unref(solver) ;
  unref(precond) ;
  unref(problem) ;
  unref(exact) ;
  unref(error) ;
  delete function1 ;
  delete function2 ;
  delete postProcessor ;
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;
  
  return 0 ;
}
