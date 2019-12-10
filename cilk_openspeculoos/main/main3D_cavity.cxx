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

// main3Dliddriven.cxx
//
// Solves the lid-driven 3D square cavity for the Navier-Stokes equations.
//
//     dv/dt + v.grad v = - grad p + 1/Re Laplacian v
//     - div v          = 0
//
// The data are:
//   . space domain: the unit square cavity ]0,1[ x ]0,1[ x ]0,1[;
//   . Re = 1000;
//   . Dirichlet conditions:
//      - on the velocity: at all times, v=0 on all faces, except on the lid
//        (top face), where it is prescribed by a regularized function
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
//   . initial conditions: at t=0, v=0 on the whole cavity.
//
// Karniadakis's method BDFn/EXn (n=1, 2, 3 or 4) is used. See Speculoos
// programmer's manual or Couzy's thesis, p 40.

//#define BDF1
#define BDF2
//#define BDF3

//#define DECOMPOSITION_TYPE UZAWA
//#define DECOMPOSITION_TYPE BP1
#define DECOMPOSITION_TYPE BP1_PC
//#define DECOMPOSITION_TYPE BP3

#ifdef PARALLEL
#define cout_root if(Mpi_rank==0) std::cout
#endif

#include "core/main.hxx"

int main (int argc, char* argv[]) {

    int      nelem1  =8, nelem2  =8, nelem3  =8 ;
    int      polydeg1=8, polydeg2=8, polydeg3=8, nstep=10 ;
    char     titre[20];
    real     rey=12000, dt=0.001, hauteur=1., longueur=1., largeur=1. ; 

    Mpi_Initialize(&argc,&argv) ;

#ifdef PARALLEL
    ParallelTimer *pt = new ParallelTimer() ;
#endif

    switch (argc) {
	case 9 :
	    nelem1   = atoi(argv[1]) ;
	    nelem2   = atoi(argv[2]) ;
	    nelem3   = atoi(argv[3]) ;
	    polydeg1 = atoi(argv[4]) ;
	    polydeg2 = atoi(argv[5]) ;
	    polydeg3 = atoi(argv[6]) ;
	    nstep    = atof(argv[7]) ;
	    dt       = atof(argv[8]) ;
	    break ;
	default :
	  //Error("main","8 arguments (nx,ny,nz,polydeg1,polydeg2,polydeg3,nstep,dt) are required") ;
	  cout << "Using default arguments" << endl;
	  break ;
    }

#ifdef PARALLEL
    pt->Start();
#endif
        
    if (Mpi_rank==0) {
	printf("nelem1  \t %2d, nelem2  \t %2d, nelem3  \t %2d\n",nelem1,nelem2,nelem3) ;
	printf("polydeg1\t %2d, polydeg2\t %2d, polydeg3\t %2d\n",polydeg1,polydeg2,polydeg3) ;
	printf("nstep\t %d\n",nstep) ;
	printf("rey\t %6f\ndt\t %f\n",rey,dt) ;
    }
    


    // ==========================
    // *** I. Mesh generation ***
    // ==========================
    Point        *p1, *p2, *p3, *p4 ;
    Point        *p5, *p6, *p7, *p8 ;
    Vertex       *v1, *v2, *v3, *v4 ;
    Vertex       *v5, *v6, *v7, *v8 ;
    Edge         *edge1, *edge2, *edge3, *edge4 ;
    Edge         *edge5, *edge6, *edge7, *edge8 ;
    Edge         *edge9, *edge10, *edge11, *edge12 ;
    Quad         *face1, *face2,*face3,*face4,*face5,*face6 ; 
    Volume       *volume ;
    Distribution *distrib1, *distrib2,*distrib3 ;
    Mesh3D       *mesh ;
    
    p1 = new Point(0.      ,0.     ,0.) ;
    p2 = new Point(longueur,0.     ,0.) ;
    p3 = new Point(longueur,hauteur,0.) ;
    p4 = new Point(0.      ,hauteur,0.) ;
    
    p5 = new Point(0.      ,0.     ,largeur) ; 
    p6 = new Point(longueur,0.     ,largeur) ; 
    p7 = new Point(longueur,hauteur,largeur) ; 
    p8 = new Point(0.      ,hauteur,largeur) ; 

    v1 = new Vertex(p1) ;
    v2 = new Vertex(p2) ;
    v3 = new Vertex(p3) ;
    v4 = new Vertex(p4) ;
    
    v5 = new Vertex(p5) ;
    v6 = new Vertex(p6) ;
    v7 = new Vertex(p7) ;
    v8 = new Vertex(p8) ;

    edge1 = new Line(v1,v2) ;   
    edge2 = new Line(v2,v3) ;   
    edge3 = new Line(v3,v4) ;   
    edge4 = new Line(v4,v1) ;   

    edge5 = new Line(v5,v6) ;   
    edge6 = new Line(v6,v7) ;   
    edge7 = new Line(v7,v8) ;   
    edge8 = new Line(v8,v5) ;   
    
    edge9  = new Line(v1,v5) ;  
    edge10 = new Line(v2,v6) ;  
    edge11 = new Line(v3,v7) ;  
    edge12 = new Line(v4,v8) ;  
    
    //--------------------------faces----------------------------
    
    face1  = new Quad(edge4,edge8,edge12,edge9 ) ;
    face2  = new Quad(edge2,edge6,edge10,edge11) ;
    face3  = new Quad(edge1,edge5,edge9 ,edge10) ;
    face4  = new Quad(edge3,edge7,edge11,edge12) ;
    face5  = new Quad(edge1,edge3,edge2 ,edge4 ) ;
    face6  = new Quad(edge5,edge7,edge8 ,edge6 ) ;
    
    //--------------------------volume----------------------------
    
    volume = new Hexa(face1,face2,face3,face4,face5,face6) ;

#ifdef PARALLEL
    cout_root << "Mesh init " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    distrib1 = new UniformDistribution(nelem1) ;
    distrib2 = new UniformDistribution(nelem2) ;
    distrib3 = new UniformDistribution(nelem3) ;

#ifdef PARALLEL
    cout_root << "Distribution " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    volume->SetDistribution(1,distrib1) ;
    volume->SetDistribution(2,distrib2) ;
    volume->SetDistribution(3,distrib3) ;

#ifdef PARALLEL
    cout_root << "volume->SetDistribution() " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    mesh = volume->GenerateMesh() ; 

#ifdef PARALLEL
    cout_root << "volume->GenerateMesh() " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    mesh->DispatchElements() ; 

#ifdef PARALLEL
    cout_root << "mesh->DispatchElements() " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    // =============================
    // *** II. Field definitions ***
    // =============================
    
    ParentElement          *parent1, *parent2, *parent3 ;
    FlatField              *coord, *pressure ;
    MortaredField          *velocity ;
    
    //    II.1 Parent element definitions and coordinate field
    //    ----------------------------------------------------

    parent1 = ParentElementGLL::GetParentElementOfDegree(polydeg1) ;
    parent2 = ParentElementGLL::GetParentElementOfDegree(polydeg2) ;
    parent3 = ParentElementGLL::GetParentElementOfDegree(polydeg3) ;
    coord   = mesh->GetCoordinates() ; 

#ifdef PARALLEL
    cout_root << "mesh->GetCoorrdinates() " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    // set polynomial degree to 'polydeg1' and 'polydeg2' in all elements
    coord->SetParentElements(parent1,1) ;
    coord->SetParentElements(parent2,2) ;
    coord->SetParentElements(parent3,3) ;

#ifdef PARALLEL
    cout_root << "coord->SetParentElements() " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    coord->SetToCoordinates() ;
#ifdef PARALLEL
    cout_root << "coord->SetToCoordinates() " << pt->GetIntermediateTime() << " sec\n" ;
#endif
  
    //    ------------------------------------------------------
    //    Velocity field (GLL identical to coordinate field)
    //    -------------------------------------------------------
    
    velocity = new MortaredField(3,coord) ;

#ifdef PARALLEL
    cout_root << "velocity = new MortaredField(3,coord) " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    //    -------------------------------------------------------
    //    II.3 Pressure field (GL)
    //    -------------------------------------------------------
    
    Vector<ParentElement*> *parents ;
    int                     polydegP, shift ;
    shift          = 2 ;  // Choice of a staggered grid: P(N)_P(N-2)
    parents        = new Vector<ParentElement*>(3) ;
    polydegP       = max(polydeg1-shift,0) ;
    parents->At(1) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
    polydegP       = max(polydeg2-shift,0) ;
    parents->At(2) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
    polydegP       = max(polydeg3-shift,0) ;
    parents->At(3) = ParentElementGL::GetParentElementOfDegree(polydegP) ;
    
    pressure = new FlatField(1,mesh,parents) ;

#ifdef PARALLEL
    cout_root << "New Pressure field (GL) " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    // =========================
    // *** III. Forcing term ***
    // =========================
    
    FlatField *forcing ;
    Function  *forcingFunction ;
    
    forcing         = NULL ; 
    forcingFunction = NULL ; 

    // =================================================
    // *** VI. Boundary conditions ***
    // =================================================
    
    Mesh2D             *boundary, *lid, *walls ;
    FlatField          *prescribedVxWalls=NULL, *prescribedVxLid ;
    FlatField          *prescribedVz, *prescribedVy ;
    Function           *functionLid ;
    DirichletCondition *dcVxWalls, *dcVxLid, *dcVy, *dcVz ;
    
    boundary = (Mesh2D*) mesh->GetSkin() ;
    lid      = (Mesh2D*) boundary->Extract(2,hauteur);
    walls    = (Mesh2D*) boundary->Except(lid);

#ifdef PARALLEL
    cout_root << "Boundary conditions " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    // VI.1. Boundary conditions on the velocity
    //------------------------------------------
    // vX=0 on the walls, vX=regularized on the lid
    // vY=0 on the whole boundary (walls + lid)
    // vZ=0 on the whole boundary (walls + lid)
    
    prescribedVxWalls = walls   ->GetCoordinates()->DuplicateEmpty(1) ;
    prescribedVxLid   = lid     ->GetCoordinates()->DuplicateEmpty(1) ;
    prescribedVy      = boundary->GetCoordinates()->DuplicateEmpty(1) ;
    prescribedVz      = boundary->GetCoordinates()->DuplicateEmpty(1) ;

#ifdef PARALLEL
    cout_root << "Boundary Conditions on the velocity " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    
    // CONDITIONS AUX LIMITES Vx SUR LA LID
    
    functionLid = new RegularizedProfile(1,18,18) ;
 
    prescribedVxLid->SetTo(functionLid) ;

    dcVxWalls = new DirichletCondition(prescribedVxWalls,1) ;
    dcVxLid   = new DirichletCondition(prescribedVxLid  ,1) ;
    dcVy      = new DirichletCondition(prescribedVy     ,2) ;
    dcVz      = new DirichletCondition(prescribedVz     ,3) ;

#ifdef PARALLEL
    cout_root << "Boundary conditions VX on the LID " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    // =========================
    // *** Integration rules ***
    // =========================
    
    FlatField *ruleMomentum, *ruleMass ;
  
    ruleMomentum = velocity->GetMain()->DuplicateEmpty(1) ;
    ruleMass     = pressure->DuplicateEmpty() ;

#ifdef PARALLEL
    cout_root << "Integration rules " << pt->GetIntermediateTime() << " sec\n" ;
#endif
    // =========================================================
    // *** V. Linear solvers (overrides the default solvers) ***
    // =========================================================
    Solver *solverH, *solverP, *solverU ;
    
    solverH = new ConjugateGradient(6000,1e-12) ;     // Helmholtz
    solverP = new ConjugateGradient(6000,1e-6) ;     // pressure (external loop)
    solverU = new ConjugateGradient(600,1e-14) ;      // J=H-1 (case Uzawa only) 
  
    solverH->SetName("Helmholtz") ;
    solverP->SetName("main pressure") ;
    solverU->SetName("internal Uzawa") ;
    
    solverH->SetVerbosity(1) ;
    solverP->SetVerbosity(1) ;
    solverU->SetVerbosity(1) ;

#ifdef PARALLEL
    cout_root << "Linear solvers " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    // ============================
    // *** VI. Time-integration ***
    // ============================
    
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
# endif

    diffusionScheme->SetTimeIncrement(dt) ;

#ifdef PARALLEL
    cout_root << "Time-integration " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    // ==================================
    // *** VII. Navier-Stokes problem ***
    // ==================================

    NavierStokes *navierStokes ;
    SteadyStokes *steadyStokes ;
    
    //    VII.1 Main settings
    //    -------------------
    
    navierStokes = new NavierStokes(mesh,DECOMPOSITION_TYPE) ;
    navierStokes->SetTimeIntegrationScheme(diffusionScheme) ;
    navierStokes->SetConvectionScheme(convectionScheme) ;
    navierStokes->SetVelocity(velocity) ;
    navierStokes->SetPressure(pressure) ;
    navierStokes->SetReynolds(rey) ;
    navierStokes->SetForcing(forcing) ;
    navierStokes->SetForcingFunction(forcingFunction) ;
    navierStokes->AddDirichletConditionV(dcVxWalls) ;
    navierStokes->AddDirichletConditionV(dcVxLid) ;
    navierStokes->AddDirichletConditionV(dcVy) ;
    navierStokes->AddDirichletConditionV(dcVz) ;
    navierStokes->SetIntegrationRuleMass(ruleMass) ;
    navierStokes->SetIntegrationRuleMomentum(ruleMomentum) ;

    steadyStokes = navierStokes->GetSteadyStokes() ;

#ifdef PARALLEL
    cout_root << "NS Problem main settings " << pt->GetIntermediateTime() << " sec\n" ;
    
#endif

    //    VII.2 Preconditioners of the 3 linear solvers
    //    ---------------------------------------------
  
    // assign the linear solvers
    steadyStokes->GetHelmholtz()->SetSolver(solverH) ;
    steadyStokes->SetSolver(solverP) ;
    if (steadyStokes->GetInternalHelmholtz() != NULL)            // case Uzawa
	steadyStokes->GetInternalHelmholtz()->SetSolver(solverU) ;
    
    Preconditioner *precondP, *precondH ;
    precondH = new DiagonalPreconditioner() ;
    precondP = new CouzyPreconditioner() ;
  
    steadyStokes->SetPreconditioner(precondP) ;
    steadyStokes->GetHelmholtz()->SetPreconditioner(precondH);
  
#ifdef PARALLEL
    cout_root << "Preconditioners of the 3 linear solvers " << pt->GetIntermediateTime() << " sec\n" ;
#endif

    //    VII.3 Initial values (and starting values, if multistep method)
    //    ---------------------------------------------------------------

    FlatField    *startV0 ;
    
    startV0 = velocity->GetMain()->DuplicateEmpty() ;           
  
    diffusionScheme->SetValues(0,startV0) ;    // v(0)=0	
  
    // ==============================
    // *** VIII. Solution process ***
    // ==============================
  
#ifdef PARALLEL
    cout_root << "Initialization time = " << pt->GetIntermediateTime() << " sec\n" ;
    pt->Stop();
#endif



    int         i ;
    double time_tmp = 0.0;
#ifdef PARALLEL
    pt->WriteStartLine();
    for (i=1 ; i<=nstep ; i++) {
	// printf("Iteration step %d\n",i) ;
	pt->Start();
	cout_root << "Loop " << i <<" Timer started \n";

//	_mpimon_start_counters();

	navierStokes->Solve(1) ;      // solve 1 step at a time
	cout_root << "--- > NS : \t" << pt->GetIntermediateTime() << " sec\n" ;
	cout_root << "--- > VF : \t" << pt->GetIntermediateTime() << " sec\n" ;
	cout_root << "Loop " << i <<" = " << pt->GetIntermediateTime() << " sec\n" ;
	cout_root << "----------------------------------------------\n";
//	_mpimon_stop_counters();
//	_mpimon_read_counters();
	time_tmp = pt->GetIntermediateTime();
	pt->WriteTimeToFile(time_tmp);
	pt->Stop();
      
    } // fin de la boucle en temps
    pt->WriteEndLine();
#else
      for (i=1 ; i<=nstep ; i++) {
	Printf("\n\nStart step %d (time = %.4f)\n\n",
	       i,diffusionScheme->GetTime()+dt) ;
	navierStokes->Solve(1) ;      // solve 1 step at a time
    }
#endif 
    // ===================
    // *** X. Cleaning ***
    // ===================
    unref(p1) ; unref(p5) ; unref(p6) ; 
    unref(p2) ; unref(p7) ; unref(p8) ;
    unref(p3) ; 
    unref(p4) ; 
    unref(v1) ;unref(v5) ;unref(v7) ;
    unref(v2) ;unref(v6) ;unref(v8) ;
    unref(v3) ;
    unref(v4) ;
    unref(edge1) ;unref(edge5) ;unref(edge8) ;unref(edge11) ;
    unref(edge2) ;unref(edge6) ;unref(edge9) ;unref(edge12) ;
    unref(edge3) ;unref(edge7) ;unref(edge10) ;
    unref(edge4) ;
    unref(face1) ; unref(face2) ; unref(face3) ;
    unref(face5) ; unref(face4) ; unref(face6) ;
    unref(distrib1) ;
    unref(distrib2) ;
    unref(distrib3) ;
    unref(mesh) ;
    unref(velocity) ;
    unref(pressure) ;
    unref(forcing) ;
    unref(forcingFunction) ;
    unref(boundary) ;
    unref(lid) ;
    unref(walls) ;
    unref(prescribedVy) ;
    unref(dcVy) ;
    unref(prescribedVz) ;
    unref(dcVz) ;
    unref(prescribedVxWalls) ;
    unref(prescribedVxLid) ;
    unref(solverH) ;
    unref(precondP) ;
    unref(solverP) ;
    unref(solverU) ;
    unref(diffusionScheme) ;
    unref(convectionScheme) ;
    unref(navierStokes) ;
    unref(startV0) ;
    ParentElement::DeleteAllInstances() ;

    Mpi_Finalize() ;

    return 0 ;
}


    
