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

// main1D_stokes.cxx
//
//
// Solves the following 1D steady Stokes problem:
//
//     lambda v = - grad p + nu Laplacian v + f
//     div v    = 0
//
// with lambda=1, nu=1 and f(x) = pi cos(pi x), on the domain ]-1,1[.
//
// The Dirichlet conditions are: v(-1) = 0 and v(1) = 0.
//
// The exact solution is: v_exact(x) = 0
//                        p_exact(x) = sin(pi x).
//
// Remarks: . the only divergence-free field v in the 1D case is v = constant.
//          . the Dirichlet conditions can be removed (the operator H is 
//            positive definite anyway).


// The following options are available:
//
// - define CHECK_FIELDS if you want to test if the solution (v,p) satisfies
//   the momentum and mass equations.
//
// - define EXPONENTIAL if you want f(x) to be replaced with f(x) =
//   2 (1-x-x*x) e**(2*x), and thus the solution p(x) by p(x) = 
//   (1-x*x) e**(2*x).
//   (and still v(x) = 0).
//

//#define CHECK_FIELDS
//#define EXPONENTIAL


#include "core/main.hxx"


int main (int argc, char* argv[])
{ // The (optional) arguments are:
  //   1. nelem:   number of elements.
  //               default = 1.
  //   2. polydeg: polynomial degree for the velocity field, in all elements.
  //               default = 2.
  //   3. shift:   degree(pressure) = degree(velocity) - shift.
  //               default = 2 (i.e., PN/PN-2 formulation).

  int nelem, polydeg, shift ;
  
  Mpi_Initialize(&argc,&argv) ;

  switch (argc) {
    case 1 :
      nelem   = 1 ;
      polydeg = 2 ;
      shift   = 2 ;
      break ;
    case 2 :
      nelem   = atoi(argv[1]) ;
      polydeg = 2 ;
      shift   = 2 ;
      break ;
    case 3 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      shift = 2 ;
      break ;
    case 4 :
      nelem   = atoi(argv[1]) ;
      polydeg = atoi(argv[2]) ;
      shift   = atoi(argv[3]) ;
      break ;
    default :
      printf("%s%s%s%s%s", "Error: valid number of arguments are:\n",
             "  0\n",
             "  1: nelem\n",
             "  2: nelem,polydeg\n",
             "  3: nelem,polydeg,shift\n") ;
      Exit(0) ;
      break ;
  }


  // =======================
  // *** Mesh generation ***
  // =======================

  Point        *p1, *p2 ;
  Vertex       *v1, *v2 ;
  Edge         *edge ;
  Distribution *distrib ;
  Mesh1D       *mesh ;

  p1   = new Point(-1.) ;
  p2   = new Point( 1.) ;

  v1   = new Vertex(p1) ;
  v2   = new Vertex(p2) ;

  edge = new Line(v1,v2) ;

  distrib = new UniformDistribution(nelem) ;
  edge->SetDistribution(1,distrib) ;

  mesh = edge->GenerateMesh() ;
  mesh->DispatchElements() ;
  mesh->PrintBarycenters() ;


  // =========================
  // *** Field definitions ***
  // =========================

  ParentElement *parent ;
  FlatField     *coord, *pressure ;
  MortaredField *velocity ;
  int           polydegP ;

  //    1. Parent element definitions and coordinate field
  //    --------------------------------------------------

  parent = ParentElementGLL::GetParentElementOfDegree(polydeg) ;
  mesh->SetDefaultParentElement(parent) ;

  coord  = mesh->GetCoordinates() ;
  printf("\nCoordinates: ") ; coord->Print() ;

  //    2. Velocity field (GLL identical to coordinate field)
  //    -----------------------------------------------------

  velocity = new MortaredField(1,coord) ;
  
  //    3. Pressure field (GL)
  //    ----------------------

  polydegP = max(polydeg-shift,0) ;
  parent   = ParentElementGL::GetParentElementOfDegree(polydegP) ;

  pressure = new FlatField(1,mesh,parent) ;


  // ============================================
  // *** Forcing term and boundary conditions ***
  // ============================================

  Mesh               *boundary ;
  FlatField          *f, *prescribedV ;
  Function           *functionForcing, *functionV, *functionP ;
  DirichletCondition *dcV ; 

  //    1. Forcing term (pi cos(pi x))
  //    ------------------------------
  
  functionForcing = new Cos(PI,PI) ;
# ifdef EXPONENTIAL
    functionForcing = new ExponentialSpecial(1) ;
# endif

  f = velocity->GetMain()->DuplicateEmpty(1) ;
  f->SetTo(functionForcing) ;
  //  printf ("\nForcing term: "); f->Print() ;

  delete functionForcing ;
  
  //    2. Boundary conditions on the velocity (set to the exact solution)
  //    ------------------------------------------------------------------

  boundary = mesh->GetSkin() ;
  printf ("\nBoundary: ") ; boundary->Print() ;

  prescribedV = new FlatField(1,coord,boundary) ;

  functionV   = new Constant(0.) ;

  functionP   = new Sin(1.,PI) ;
# ifdef EXPONENTIAL
    functionP = new ExponentialSpecial(0) ;
# endif

  prescribedV->SetTo(functionV) ;
  //  printf ("\nPrescribed on v: "); prescribedV->Print() ;

  dcV = new DirichletCondition(prescribedV,1) ;


  // =========================
  // *** Integration rules ***
  // =========================

  FlatField *ruleMomentum, *ruleMass ;

  ruleMomentum = velocity->GetMain()->DuplicateEmpty(1) ;
  ruleMass     = pressure->DuplicateEmpty() ;


  // ==============
  // *** Solver ***
  // ==============

  // only the solver for the pressure equation. The solver for the momentum 
  // equation and the solver for the inner loop of the Uzawa operator are
  // the default ones (and given in class SteadyStokes).

  Solver *solver ;

  //  solver = new DirectSolver(true) ;
  solver = new ConjugateGradient(300,1e-12) ;
  solver->SetVerbosity(2) ;


  // ===============
  // *** Problem ***
  // ===============

  SteadyStokes *problem ;
  real         lambda, nu ;

  lambda = 1. ;
  nu     = 1. ;

  problem = new SteadyStokes(mesh,UZAWA) ;
  problem->SetLambda(lambda) ;
  problem->SetNu(nu) ;
  problem->SetVelocity(velocity) ;
  problem->SetPressure(pressure) ;
  problem->SetIntegrationRuleMomentum(ruleMomentum) ;
  problem->SetIntegrationRuleMass(ruleMass) ;
  problem->SetForcing(f) ;
  problem->AddDirichletConditionV(dcV) ;
  problem->SetSolver(solver) ;


  // ========================
  // *** Solution process ***
  // ========================

  printf("\nStart solving\n") ;

  problem->Solve() ;

  //  printf("\nSolution: velocity: ") ; problem->GetVelocity()->Print() ;
  //  printf("\nSolution: pressure: ") ; problem->GetPressure()->Print() ;


  // ======================
  // *** Postprocessing ***
  // ======================

  FlatField *vel, *errorV, *errorP ;

  vel    = velocity->GetMain() ;

  errorV = vel->DuplicateEmpty() ;
  errorP = pressure->DuplicateEmpty() ;

  errorV->SetTo(functionV) ;
  errorP->SetTo(functionP) ;
  errorP->SetToZeroMean(1) ;
  //  printf("\nExact pressure: "); errorP->Print() ;

  errorV->Subtract(vel) ;                           // errorV = v(exact) - v
  errorP->Subtract(pressure) ;                      // errorP = p(exact) - p

# ifdef CHECK_FIELDS

  MortaredField *gradP, *vwork ;
  FlatField     *v, *Hv, *Lv, *Gp, *Bf, *Dv, *errorMomentum, *errorUzawa,
                *ruleMoment, *ruleMass, *lhs, *rhs ;

  v             = velocity->GetMain() ;
  Hv            = v->DuplicateEmpty() ;
  Lv            = v->DuplicateEmpty() ;
  Gp            = v->DuplicateEmpty() ;
  Bf            = v->DuplicateEmpty() ;
  errorMomentum = v->DuplicateEmpty() ;
  
  vwork         = velocity->DuplicateEmpty() ;
  gradP         = velocity->DuplicateEmpty() ;
  lhs           = pressure->DuplicateEmpty() ;
  rhs           = pressure->DuplicateEmpty() ;
  Dv            = pressure->DuplicateEmpty() ;
  errorUzawa    = pressure->DuplicateEmpty() ;
 
  ruleMoment    = problem->GetIntegrationRuleMoment() ;
  ruleMass      = problem->GetIntegrationRuleMass() ;

  // 1. Check momentum equation H v + Grad(p) - B f = 0

  Hv->SetToWeakIdentity(v,ruleMoment) ;
  Hv->Multiply(lambda) ;
  Lv->SetToWeakLaplacian(v,ruleMoment) ;
  Hv->Add(Lv,-nu) ;

  Gp->SetToWeakGradient(pressure,ruleMass) ;

  Bf->SetToWeakIdentity(f,ruleMoment) ;

  //    printf("\nH v : ")    ; Hv->Print() ;
  //    printf("\nGrad(p): ") ; Gp->Print() ;
  //    printf("\nB f : ")    ; Bf->Print() ;

  errorMomentum->CopyFrom(Hv) ;
  errorMomentum->Add(Gp) ;
  errorMomentum->Subtract(Bf) ;
  printf("\nError on momentum equation: H v + Grad(p) - B f :\n") ;
  //  errorMomentum->Print() ;

  // 2. Check Uzawa equation -D J D(T) p - D v~ = 0

  gradP->GetMain()->SetToWeakGradient(pressure,ruleMass) ;
  problem->ApplyJ(gradP,vwork) ;
  lhs->SetToWeakDivergence(vwork->GetMain(),ruleMass);
  rhs->SetToWeakDivergence(problem->GetVtilde(),ruleMass) ;

  //   printf("\nLHS of Uzawa: ") ; lhs->Print() ;
  //   printf("\nRHS of Uzawa: ") ; rhs->Print() ;

  errorUzawa->CopyFrom(lhs) ;
  errorUzawa->Subtract(rhs) ;
  printf("\nError on Uzawa: lhs - rhs :\n") ; errorUzawa->Print() ;

  // 3. Complementary check: Div v = 0

  Dv->SetToWeakDivergence(v,ruleMass) ;
  //  printf("\nDiv(v) : ") ; Dv->Print() ;

# endif

  // Velocity and pressure

  real normV, normP ;

  //  printf("\nError on v: ") ; errorV->Print() ;
  //  printf("\nError on p: ") ; errorP->Print() ;

  normV = errorV->GetH1Norm() ;
  normP = errorP->GetL2Norm() ;
  
  printf("\nH1 norm of error field on v:  %12.4E",normV) ;
  printf("\nL2 norm of error field on p:  %12.4E\n\n",normP) ;
  
  // Plots
 
  FILE *plotFile ;

  plotFile = Fopen("./stokes1Dsteady.dat","a+") ;
  fprintf(plotFile,"%d %e %e\n",polydeg,normV,normP) ;
  fclose(plotFile) ;


  // ================
  // *** Cleaning ***
  // ================

  unref(p1) ;
  unref(p2) ;
  unref(v1) ;
  unref(v2) ;
  unref(edge) ;
  unref(distrib) ;
  unref(mesh) ;
  unref(velocity) ;
  unref(pressure) ;
  unref(boundary) ;
  unref(f) ;
  unref(prescribedV) ;
  unref(dcV) ;
  unref(ruleMomentum) ;
  unref(ruleMass) ;
  unref(functionV) ;
  unref(functionP) ;
  unref(solver) ;
  unref(problem) ;
/*
  unref(errorV) ;
  unref(errorP) ;
*/
  ParentElement::DeleteAllInstances() ;

  Mpi_Finalize() ;

  return 0 ;
}
