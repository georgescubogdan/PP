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

// timinteg.cxx

#include "core/timinteg.hxx"
#include "core/problem.hxx"
#include "core/field.hxx"
#include "core/solver.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include "core/parallel.hxx"
#include <cstring>

//___________________________ TimeIntegrationScheme ___________________________


TimeIntegrationScheme :: TimeIntegrationScheme ()
    : GarbagedObject ()
{
    // Constructor. Initializes a new scheme.

    currentTimeStep       = 0 ;
    time                  = 0. ;
    dt                    = UNINITIALIZED_REAL ;
    semiDiscreteProblem   = NULL ;
    solution              = NULL ;
    pastSolutions         = NULL ;
    additionalWeakHistory = NULL ;
    starter               = NULL ;
    useStarter            = true ;
    order                 = UNINITIALIZED_INTEGER ;
}


TimeIntegrationScheme :: ~TimeIntegrationScheme ()
{
    // Destructor.

    int i ;

    unref(semiDiscreteProblem) ;
    unref(solution) ;
    unref(additionalWeakHistory) ;
    unref(starter) ;

    if (pastSolutions) {
	for (i=1 ; i<=pastSolutions->GetSize() ; i++)
	    unref(pastSolutions->At(i)) ;
	delete pastSolutions ;
    }
}


void TimeIntegrationScheme :: ComputePastSequence (FlatField*)
{
    // Implemented by subclasses.

    ImplementedBySubclasses("TimeIntegrationScheme::ComputePastSequence(y)") ;
}

int TimeIntegrationScheme :: GetOrder()
{
    return order;
}

Field* TimeIntegrationScheme :: GetSolution ()
{
    // Returns the current value of the solution field.
    // If the current time step is n, returns u(n).

    return solution ;
}

Field* TimeIntegrationScheme :: GetPastSolution (int q)
{
    // Returns the current value of the solution field.
    // If the current time step is n, returns u(n).

#ifdef REQUIRE
    Require("order greater or equal to q", q <= pastSolutions->GetSize());
    InMethod("TimeIntegrationScheme :: GetPastSolution (int q)");
#endif

    if (pastSolutions) 
	return pastSolutions->At(q) ;
    else
	return NULL;
}

real TimeIntegrationScheme :: GetTime ()
{
    // Returns the current date.

    return time ;
}
real TimeIntegrationScheme :: GetTimeIncrement ()
{
    // Returns the current date.

    return dt ;
}

void TimeIntegrationScheme :: PrintInProc0Past(const char *nom)
{ 

    char liste[7];
    int n, i;
    char *nomx, *nomy, *nomz;
    int nbcomp ;

    nbcomp = pastSolutions->At(1)->GetNbComponents();

    nomx = new char[strlen(nom)+4];
    nomy = new char[strlen(nom)+4];
    if (nbcomp==3)  
	nomz = new char[strlen(nom)+4];
  
    n = pastSolutions->GetSize() ;
  
    strcpy(liste,"012345");
    strcpy(nomx,"ux");
    strcpy(nomy,"uy");
    if (nbcomp==3) 
	strcpy(nomz,"uz");
  
    strcat(nomx,nom);
    strcat(nomy,nom);
    if (nbcomp==3) 
	strcat(nomz,nom);
  
    for (i=1 ; i<=n ; i++)  
    {
	nomx[6] = liste[i];
	nomy[6] = liste[i];
	nomx[7] = '\0';
	nomy[7] = '\0';
	pastSolutions->At(i)->PrintInProc0(nomx,1);
	pastSolutions->At(i)->PrintInProc0(nomy,2);
      
	if (nbcomp==3){
	    nomz[6] = liste[i];
	    nomz[7] = '\0';
	    pastSolutions->At(i)->PrintInProc0(nomz,3);
	}
    }
  
    delete(nomx);
    delete(nomy);
    if (nbcomp==3)  
	delete(nomz);
}

void TimeIntegrationScheme :: ReadFromProc0Past(const char *nom)
{
  
    char liste[7];
    char *nomx, *nomy, *nomz;
    int n, i, nbcomp;
  
    nbcomp = pastSolutions->At(1)->GetNbComponents();
  
    nomx = new char[strlen(nom)+5];
    nomy = new char[strlen(nom)+5];
    if (nbcomp==3) 
	nomz = new char[strlen(nom)+5];
  
    n = pastSolutions->GetSize() ;
  
    strcpy(liste,"012345");
    strcpy(nomx,"ux");
    strcpy(nomy,"uy");
    if (nbcomp==3) 
	strcpy(nomz,"uz");
  
    strcat(nomx,nom);
    strcat(nomy,nom);
    if (nbcomp==3) 
	strcat(nomz,nom);
  
    for (i=1 ; i<=n ; i++)  
    {
	nomx[6] = liste[i];
	nomy[6] = liste[i];
      
	nomx[7] = '\0';
	nomy[7] = '\0';
      
	pastSolutions->At(i)->ReadFromProc0(nomx,1);
	pastSolutions->At(i)->ReadFromProc0(nomy,2);
      
	if (nbcomp==3){
	    nomz[6] = liste[i];
//	nomz[7] = '0';
	    nomz[7] = '\0';
	    nomz[8] = '\0';
	    pastSolutions->At(i)->ReadFromProc0(nomz,3);
	}
    }
  
    delete(nomx);
    delete(nomy);
    if (nbcomp==3)  
	delete(nomz);
}

void TimeIntegrationScheme :: SetAdditionalWeakHistory (FlatField* hist)
{
    // At every step the receiver will add 'hist' to its weak history field.
    // 'hist' may be NULL.

    ref(hist) ;
    unref(additionalWeakHistory) ;
    additionalWeakHistory = hist ;
}


void TimeIntegrationScheme :: SetCurrentTimeStep (int n)
{
    // Sets to 'n' the number of the current time step.

    currentTimeStep = n ;
}


void TimeIntegrationScheme :: SetField (Field* x)
{
    // Sets the solution field of the receiver to 'x'.
    // The solution field of the semi-discrete problem and the current solution
    // of the receiver are both 'x'.
    // The values of 'x' may contain either
    // - adequate initial values, or
    // - garbage, in which case they can be subsequently initialized by invoking
    //   method SetValues(0,x0).
    // The latter case is recommended, because it applies also to initializing
    // the field at starting steps (for a non self-starting scheme).

    int i ;

    ref(x) ;
    unref(solution) ;
    solution = x ;

    if (semiDiscreteProblem)
	semiDiscreteProblem->SetTemporalUnknown(x) ;

    // allocate the previous solutions
    for (i=1 ; i<=pastSolutions->GetSize() ; i++) {
	unref(pastSolutions->At(i)) ;
	pastSolutions->At(i) = x->GetMain()->DuplicateEmpty() ;
    }

    currentTimeStep = 0 ;
}
void TimeIntegrationScheme ::PrintPastSequence ()
{
    // default does nothing
}

void TimeIntegrationScheme :: SetSemiDiscreteProblem(SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the semi-discrete problem that the receiver will solve at
    // each time step.

# ifdef REQUIRE
    Require("has no semi-discrete problem yet", semiDiscreteProblem == NULL) ;
    Require("same mesh", solution == NULL || 
	    solution->GetMesh() == prob->GetMesh()) ;
    InMethod("TimeIntegrationScheme::SetSemiDiscreteProblem(prob)") ;
# endif

    ref(prob) ;
    unref(semiDiscreteProblem) ;
    semiDiscreteProblem = prob ;

    if (solution)
	prob->SetTemporalUnknown(solution) ;
}


void TimeIntegrationScheme :: SetValues (int i, FlatField* x)
{
    // Sets to 'x' the past values of the unknown field u at step n-'i'.
    // For example, i=0,1,2 initializes respectively u(n), u(n-1) and u(n-2) for
    // next time step n+1.
    // For example, in order to give the initial values u(0) for a self-starting
    // scheme such as BDF1 or Crank-Nicolson, invoke this method with 'i'=0.

# ifdef REQUIRE
    Require("valid argument 'i'", i >= 0 && i < pastSolutions->GetSize()) ;
    Require("type of solution field already set", solution != NULL) ;
    Require("same mesh", solution==NULL || x->GetMesh()==solution->GetMesh()) ;
    Require("same interpolation", solution==NULL ||
	    x->HasSameInterpolationAs(solution->GetMain()));
    InMethod("TimeIntegrationScheme::SetValues(i,x)") ;
# endif

    if (i == 0)
	solution->GetMain()->CopyFrom(x) ;
    else
	pastSolutions->At(i)->CopyFrom(x) ;

    // note. updating will take place before the computations of next time step
    //       begin
}


void TimeIntegrationScheme :: SetTime (real t)
{
    // Sets the current time to 't'.

    time = t ;
}


void TimeIntegrationScheme :: SetTimeIncrement (real deltaT)
{
    // Sets the time increment of the receiver to 'deltaT'.

# ifdef REQUIRE
    Require("valid argument 'deltaT'", deltaT > 0.) ;
    InMethod("TimeIntegrationScheme::SetTimeIncrement(deltaT)") ;
# endif

    dt = deltaT ;
}


void TimeIntegrationScheme :: UseStarter (boolean b)
{
    // If 'b' is True, the starting time step(s) will be performed by the
    // default starting scheme (defined in method GetStarter()).
    // If 'b' is False, the starting time step(s) will be performed using the
    // values which must be supplied by means of method SetValues.
    // Ignored for a self-starting scheme (e.g., BDF1, Crank-Nicolson).

    useStarter = b ;
}


//______________________________ AdamsBashforth _______________________________


AdamsBashforth :: AdamsBashforth (int ord)
    : TimeIntegrationScheme()
{
    // Constructor. Initializes the receiver to an Adams-Bashforth scheme of
    // order 'ord'.

# ifdef REQUIRE
    Require("valid order", ord > 0 && ord <= 4) ;
    InMethod("AdamsBashforth::AdamsBashforth(ord)") ;
# endif

    int i ;

    order        = ord ;
    coefficients = NULL ;

    pastSolutions = new Vector<FlatField*>(1) ;
    pastSolutions->At(1) = NULL ;

    pastFunctionals = new Vector<FlatField*>(order) ;
    for (i=1 ; i<=order ; i++)
	pastFunctionals->At(i) = NULL ;

    SetCoefficients() ;

}


AdamsBashforth :: ~AdamsBashforth ()
{
    // Destructor.

    int i ;

    delete coefficients ;

    for (i=1 ; i<=pastFunctionals->GetSize() ; i++)
	unref(pastFunctionals->At(i)) ;
    delete pastFunctionals ;
}


void AdamsBashforth :: ComputePastSequence (FlatField* y)
{
    // Initializes 'y' to 
    //   Sum[coeff(i) * F(u(i))], i = n,...,n-order+1,
    // i.e., to the weak history vector at the current time step.
    // Typically used for integrating the nonlinear part of a Navier-Stokes
    // problem.

    int i ;

    y->SetValues(ZERO) ;
    for (i=1 ; i<=order ; i++)
	y->Add(pastFunctionals->At(i),coefficients->At(i)) ;
}


FlatField* AdamsBashforth :: ComputeStrongHistory ()
{
    // Returns the "strong" part of the history vector at the current time step
    // n+1:
    //   u(n) / dt .
    // In order to save space, the answer overwrites u(n).

    pastSolutions->At(1)->Multiply(ONE/dt) ;

    return pastSolutions->At(1) ;
}

  
FlatField* AdamsBashforth :: ComputeWeakHistory ()
{
    // Returns the "weak" part of the history vector at the current time step
    // n+1:
    //   Sum[coeff(i) * F(u(i))], i = n,...,n-order+1.
    // In order to save space, the answer overwrites the oldest past 
    // functional F(u(n-order+1)).

    FlatField *answer ;
    int       i ;

    answer = pastFunctionals->At(order) ;
    answer->Print() ;
    cout << "******************* pastfunctionals1 (order) ***********" << endl ;
    answer->Multiply(coefficients->At(order)) ;

    for (i=1 ; i<order ; i++)
	answer->Add(pastFunctionals->At(i),coefficients->At(i)) ;
    answer->Print() ;
    cout << "******************* pastfunctionals2 (order) ***********" << endl ;
    return answer ;
}
  
void AdamsBashforth :: DoOneStep ()
{
    // Updates the temporal unknown by one time step, i.e., performs
    //   u(n+1)/dt = u(n)/dt + Sum[coeff(i) * F(u(i))], i = n,...,n-order+1.

# ifdef REQUIRE
    Require("has a semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has initial values", solution != NULL) ;
    Require("has dt", dt != UNINITIALIZED_REAL) ;
    InMethod("AdamsBashforth::DoOneStep()") ;
# endif

    FlatField *weak, *strong ;
 
    if (currentTimeStep == 0)                   // perform step 0! 
	semiDiscreteProblem->GetFunctional(solution->GetMain(),time, 
					   pastFunctionals->At(order)) ; 
    UpdateTime() ;
 
    if (currentTimeStep < order && useStarter){
	// starting steps: 1,...,order-1
	GetStarter() ;
	starter->DoOneStep() ;
	semiDiscreteProblem->GetFunctional(solution->GetMain(),time,
					   pastFunctionals->At(order)) ;
	if (currentTimeStep == order-1) {
	    unref(starter) ;
	    starter = NULL ;
	    semiDiscreteProblem->SetToExplicit() ;
	}
    }
    else {                                      // normal step 
	// set lambda
	semiDiscreteProblem->SetLambda(ONE/dt) ;

	// compute the history field Hist(n+1)
	weak   = ComputeWeakHistory() ;                   // overwrites oldest F(i)
	weak->Print() ;
	cout << "******************** WEAK ************************" << endl ;
	strong = ComputeStrongHistory() ;                 // overwrites u(n)
	strong->Print() ;
	cout << "******************** STRONG ************************" << endl ;
	semiDiscreteProblem->SetWeakHistory(weak) ;
	semiDiscreteProblem->SetStrongHistory(strong) ;

	// get u(n+1) ;
	semiDiscreteProblem->Solve() ;                    // explicit solving
	solution->GetMain()->Print() ;
	cout << "******************** Solution ************************" << endl ;

	// compute F(u(n+1)) (overwrites oldest functional)
	semiDiscreteProblem->GetFunctional(solution->GetMain(),time,
					   pastFunctionals->At(order)) ;
    }
}


TimeIntegrationScheme* AdamsBashforth :: GetStarter ()
{
    // Returns the starting scheme of the receiver. Creates it if needed (i.e.,
    // if order > 1) and if the user has not provided one.

# ifdef REQUIRE
    Require("order > 1", order > 1) ;
    Require("has initial values", solution != NULL) ;
    Require("has semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has time increment dt", dt != UNINITIALIZED_REAL) ;
    InMethod("AdamsBashforth::GetStarter()") ;
# endif

    if (starter == NULL) {
	if (order == 2)
	    starter = new ThetaMethod(0.5) ;
	else
	    starter = new RungeKutta(4) ;
	starter->SetSemiDiscreteProblem(semiDiscreteProblem) ;
	starter->SetField(solution) ;
	starter->SetTimeIncrement(dt) ;
    }

    return starter ;
}


void AdamsBashforth :: SetCoefficients ()
{
    // Initializes the coefficients of the receiver.

    if (order == 1) {
	coefficients = new RealVector(1) ;
	coefficients->At(1) = 1. ;
    }
    else if (order == 2) {
	coefficients = new RealVector(2) ;
	coefficients->At(1) =  1.5 ;
	coefficients->At(2) = -0.5 ;
    }
    else if (order == 3) {
	coefficients = new RealVector(3) ;
	coefficients->At(1) =  23. / 12. ;
	coefficients->At(2) = -16. / 12. ;
	coefficients->At(3) =   5. / 12. ;
    }
    else if (order == 4) {
	coefficients = new RealVector(4) ;
	coefficients->At(1) =  55. / 24. ;
	coefficients->At(2) = -59. / 24. ;
	coefficients->At(3) =  37. / 24. ;
	coefficients->At(4) = - 9. / 24. ;
    }
    else
	Error("AdamsBashforth::SetCoefficients","unexpected order") ;
}


void AdamsBashforth :: SetField (Field* x)
{
    // Sets the solution field of the receiver to 'x'.
    // See the description in version of class TimeIntegrationScheme.

    int i ;
  
    TimeIntegrationScheme::SetField(x) ;

    // allocate the past functionals
    for (i=1 ; i<=pastFunctionals->GetSize() ; i++) {
	unref(pastFunctionals->At(i)) ; 
	pastFunctionals->At(i) = x->GetMain()->DuplicateEmpty() ;
    }
}


void AdamsBashforth :: SetSemiDiscreteProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the semi-discrete problem that the receiver will solve at
    // each time step.

# ifdef REQUIRE
    Require("has no semi-discrete problem yet", semiDiscreteProblem == NULL) ;
    Require("same mesh", solution == NULL || 
	    solution->GetMesh() == prob->GetMesh()) ;
    InMethod("AdamsBashforth::SetSemiDiscreteProblem(prob)") ;
# endif

    TimeIntegrationScheme::SetSemiDiscreteProblem(prob) ;

    semiDiscreteProblem->SetToExplicit() ;

}


void AdamsBashforth :: SetTimeIncrement (real deltaT)
{
    // Sets (or resets) the time increment of the receiver to 'deltaT'.
    // The time increment cannot be modified for orders >= 2, because the past
    // values become inconsistent.

# ifdef REQUIRE
    Require("time increment set once", dt == UNINITIALIZED_REAL || order <= 1) ;
    Require("'deltaT' > 0", deltaT > 0.) ;
    InMethod("ThetaMethod::SetTimeIncrement(deltaT)") ;
# endif

    dt = deltaT ;
}


void AdamsBashforth :: UpdateTime ()
{
    // Advances the date by one time increment and updates the past values.

    FlatField *Fn ;
    int       i, m ;

    currentTimeStep++ ;
    time = time + dt ;

    // permute past F(i). New sequence: F(n), F(n-1), etc
    m  = pastFunctionals->GetSize() ;
    Fn = pastFunctionals->At(m) ;    // F(n) was stored there at end of last step
    for (i=m ; i>=2 ; i--)
	pastFunctionals->At(i) = pastFunctionals->At(i-1) ; 
    pastFunctionals->At(1) = Fn ;

    // update 'pastSolution'
    pastSolutions->At(1)->CopyFrom(solution->GetMain()) ;
}
 

//___________________________________ BDF _____________________________________


BDF :: BDF (int ord)
    : TimeIntegrationScheme ()
{
    // Constructor. Initializes the receiver to a backward-difference scheme of
    // order 'ord'.

# ifdef REQUIRE
    Require("valid order", ord > 0 && ord <= 4) ;
    InMethod("BDF::BDF(ord)") ;
# endif

    int i ;

    order         = ord ;
    a0            = UNINITIALIZED_REAL ;
    coefficients  = NULL ;

    pastSolutions = new Vector<FlatField*>(order) ;
    for (i=1 ; i<=order ; i++)
	pastSolutions->At(i) = NULL ;

    SetCoefficients() ;
}


BDF :: ~BDF ()
{
    // Destructor.

    delete coefficients ;
}


FlatField* BDF :: ComputeStrongHistory ()
{
    // Returns the "strong" part of the history vector at the current time step
    // n+1:
    //   (-1/dt) Sum [coeff(i) * u(i)], i = n,...,n-order+1.
    // In order to save space, the answer overwrites the oldest past solution
    // u(n-order+1).

    FlatField *answer ;
    int       i ;

    answer = pastSolutions->At(order) ;                        // overwritten
    answer->Multiply(-coefficients->At(order)/dt) ;
    for (i=1 ; i<order ; i++)
	answer->Add(pastSolutions->At(i),-coefficients->At(i)/dt) ;
  
    return answer ;
}

void BDF :: RetrievePastSolutionAtOrder()
{
    FlatField *answer ;
    int       i ;

    answer = pastSolutions->At(order) ;                        // overwritten
    for (i=1 ; i<order ; i++)
	answer->Add(pastSolutions->At(i),coefficients->At(i)/dt) ;
  
    answer->Multiply(-dt/coefficients->At(order)) ;
}

void BDF :: DoOneStep ()
{
    // Advances the temporal unknown by one time step, i.e., solves for u(n+1):
    //   a_0/dt u(n+1) = F(n+1) - 1/dt Sum {a_1 u(n) + a_2 u(n-1) + ...}.
    // A scheme of order n greater than 1 is started with a BDF scheme of order
    // n-1.

# ifdef REQUIRE
    Require("has a semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has initial values", solution != NULL) ;
    Require("has dt", dt != UNINITIALIZED_REAL) ;
    InMethod("BDF::DoOneStep()") ;
# endif

    FlatField *hist ;

    UpdateTime() ;

    if (currentTimeStep < order && useStarter) {  // starting step 1,...,order-1
	GetStarter() ;
	starter->DoOneStep() ;
	if (currentTimeStep == order-1) {
	    unref(starter) ;
	    starter = NULL ;
	    semiDiscreteProblem->SetWeakHistory(NULL) ;
	}
    }

    else {                                        // normal step 
	// set mode
	semiDiscreteProblem->SetToImplicit() ;

	// set lambda
	semiDiscreteProblem->SetLambda(a0/dt) ;

	// compute the strong history field (overwrites the oldest past solution)
	hist = ComputeStrongHistory() ;
	semiDiscreteProblem->SetStrongHistory(hist) ;

	// add the additional weak history field, if any
	semiDiscreteProblem->SetWeakHistory(additionalWeakHistory) ;

	// get u(n+1) ;
	semiDiscreteProblem->Solve() ;
    }
}


TimeIntegrationScheme* BDF :: GetStarter ()
{
    // Returns the starting scheme of the receiver. Creates it if needed (i.e.,
    // if order > 1) and if the user has not provided any.

# ifdef REQUIRE
    Require("order > 1", order > 1) ;
    Require("has initial values", solution != NULL) ;
    Require("has semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has time increment dt", dt != UNINITIALIZED_REAL) ;
    InMethod("BDF::GetStarter()") ;
# endif

    if (starter == NULL) {
	// modif E.Perchat Juillet 2001
//     if (order == 2)
//       starter = new ThetaMethod(0.5) ;
//     else if (order == 3)
//       starter = new ThetaMethod(0.5) ;
//     else
//       starter = new BDF(1) ;
	if (order == 2)
	    starter = new BDF(1);
	else if (order == 3)
	    starter = new BDF(2); // coment mettre RK 2 !!!
	//      starter = new RungeKutta(2); // coment mettre RK 2 !!!
	else
	    starter = new BDF(4) ; //BDF(3) ;// coment mettre RK 4 !!!
    
	starter->SetSemiDiscreteProblem(semiDiscreteProblem) ;
	starter->SetField(solution) ;
	starter->SetTimeIncrement(dt) ;
	starter->SetAdditionalWeakHistory(additionalWeakHistory) ;
    }

    return starter ;
}

real* BDF :: GetCoefficients ()
{
    real* coeff;
    coeff = new real[order + 1];

    coeff[0] = a0;
    for (int i=1; i<=order; i++)
	coeff[i] = coefficients->At(i);

    return coeff;
}

void BDF :: SetCoefficients ()
{
    // Initializes the coefficients of the receiver.

    if (order == 1) {
	a0           = 1. ;
	coefficients = new RealVector(1) ;
	coefficients->At(1) = -1. ;
    }
    else if (order == 2) {
	a0           = 1.5 ;
	coefficients = new RealVector(2) ;
	coefficients->At(1) = -2.  ;
	coefficients->At(2) =  0.5 ;
    }
    else if (order == 3) {
	a0           = 11./6. ;
	coefficients = new RealVector(3) ;
	coefficients->At(1) = -3.  ;
	coefficients->At(2) =  1.5 ;
	coefficients->At(3) = -1. / 3. ;
    }
    else if (order == 4) {
	a0           = 25./12. ;
	coefficients = new RealVector(4) ;
	coefficients->At(1) = -4. ;
	coefficients->At(2) =  3. ;
	coefficients->At(3) = -4. / 3. ;
	coefficients->At(4) =  0.25 ;
    }
    else
	Error("BDF::SetCoefficients","unexpected order") ;
}


void BDF :: SetTimeIncrement (real deltaT)
{
    // Sets (or resets) the time increment of the receiver to 'deltaT'.
    // The time increment cannot be modified for orders >= 2, because the past
    // values become inconsistent.

# ifdef REQUIRE
    Require("time increment set once", dt == UNINITIALIZED_REAL || order <= 1) ;
    Require("'deltaT' > 0", deltaT > 0.) ;
    InMethod("ThetaMethod::SetTimeIncrement(deltaT)") ;
# endif

    dt = deltaT ;
}


void BDF :: UpdateTime ()
{
    // Advances the date by one time increment and updates the past values.

    FlatField *temp ;
    int       i, n ;

    currentTimeStep++;
    time = time + dt ;

    semiDiscreteProblem->SetTime(time) ;

    n    = pastSolutions->GetSize() ;
    temp = pastSolutions->At(n) ;
    for (i=n ; i>=2 ; i--)                              // do permutations to
	pastSolutions->At(i) = pastSolutions->At(i-1) ;   // avoid re-allocations
    pastSolutions->At(1) = temp ;
    pastSolutions->At(1)->CopyFrom(solution->GetMain()) ;
}


//_______________________________ Extrapolation _______________________________


Extrapolation :: Extrapolation (int np)
    : TimeIntegrationScheme ()
{
    // Constructor. Initializes the receiver to an 'np'-point extrapolation 
    // scheme.

# ifdef REQUIRE
    Require("valid number of points", np >= 1 && np <= 4) ;
    InMethod("Extrapolation::Extrapolation(np)") ;
# endif

    int i ;

    order        = np;
    if (Mpi_rank==0) printf(" New Extrapolation of order %d \n", order) ;
    nbPoints     = np ;
    coefficients = NULL ;

    pastSolutions = new Vector<FlatField*>(nbPoints) ;
    for (i=1 ; i<=pastSolutions->GetSize() ; i++)
	pastSolutions->At(i) = NULL ;

    SetCoefficients() ;
}


Extrapolation :: ~Extrapolation ()
{
    // Destructor.

    delete coefficients ;
}

void Extrapolation :: ComputePastSequence (FlatField* y)
{
    // Initializes 'y' to the sum of the weighted past values.

    int i ;
  
    y->SetValues(ZERO) ;
    for (i=1 ; i<=nbPoints ; i++)
	y->Add(pastSolutions->At(i),coefficients->At(i)) ;
}

 
void Extrapolation :: DoOneStep ()
{
    // Not implemented.

    NotImplemented("Extrapolation::DoOneStep()") ;
}


void Extrapolation :: SetCoefficients ()
{
    // Initializes the coefficients of the receiver.

    if (nbPoints == 1) {
	coefficients = new RealVector(1) ;
	coefficients->At(1) = 1. ;
    }
    else if (nbPoints == 2) {
	coefficients = new RealVector(2) ;
	coefficients->At(1) =  2. ;
	coefficients->At(2) = -1. ;
    }
    else if (nbPoints == 3) {
	coefficients = new RealVector(3) ;
	coefficients->At(1) =  3. ;
	coefficients->At(2) = -3. ;
	coefficients->At(3) =  1. ;
    }
    else if (nbPoints == 4) {
	coefficients = new RealVector(4) ;
	coefficients->At(1) =  4. ;
	coefficients->At(2) = -6. ;
	coefficients->At(3) =  4. ;
	coefficients->At(4) = -1. ;
    }
    else
	Error("Extrapolation::SetCoefficients","unexpected number of points") ;
}


void Extrapolation :: SetField (Field* x)
{
    // Sets the interpolation of the past solutions of the receiver to that of
    // 'x'.
    // As opposed to other time-integration-scheme classes, no reference to 'x'
    // is kept here.

    int i ;

    for (i=1 ; i<=pastSolutions->GetSize() ; i++) {
	unref(pastSolutions->At(i)) ;
	pastSolutions->At(i) = x->GetMain()->DuplicateEmpty() ;
    }
}

void Extrapolation :: PrintPastSequence ()
{

    for (int i=1 ; i<=nbPoints ; i++)
    {
	printf(" \n Past solution numero %d  coeff = %f\n",i,coefficients->At(i));
	pastSolutions->At(i)->Print();
    }
  
}

void Extrapolation :: SetSemiDiscreteProblem (SemiDiscreteProblem*)
{
    // Not implemented.

    NotImplemented("Extrapolation::SetSemiDiscreteProblem(prob)") ;
}


void Extrapolation :: SetValues (int i, FlatField* x)
{
    // Sets to 'x' the past values of the receiver at step n-'i'.
    // For example, for the EX2 scheme, i=1 (resp. i=3) corresponds to the newest
    // (resp. oldest) fields of the sequence.
    // As opposed to other time-integration scheme classes, the caller should 
    // not use i=0!

# ifdef REQUIRE
    Require("valid argument 'i'", i >= 1 && i <= pastSolutions->GetSize()) ;
    Require("type of solution field already set", pastSolutions->At(1) != NULL) ;
    Require("same mesh", x->GetMesh() == pastSolutions->At(1)->GetMesh()) ;
    Require("same interpolation",x->HasSameInterpolationAs(
		pastSolutions->At(1))) ;
    InMethod("Extrapolation::SetValues(i,x)") ;
# endif

    if (i <= pastSolutions->GetSize())
	pastSolutions->At(i)->CopyFrom(x) ;
}

void Extrapolation :: UpdateTime ()
{
    // Updates the past solutions of the receiver. At method's exit, the newest
    // past solution u(n) contains garbage.
    // As opposed to other time-integration-scheme classes, attribute 
    // <{solution}> is not used.

    FlatField *oldest ;
    int       i ;

    // use cyclic permutations for economy

    oldest = pastSolutions->At(nbPoints) ;
    for (i=nbPoints ; i>=2 ; i--)
	pastSolutions->At(i) = pastSolutions->At(i-1) ;
    pastSolutions->At(1) = oldest ;
}


//________________________________ RungeKutta _________________________________


RungeKutta :: RungeKutta (int ord)
    : TimeIntegrationScheme ()
{
    // Constructor. Initializes the receiver to a Runge-Kutta scheme of order 
    // 'ord'.

# ifdef REQUIRE
    Require("valid order", ord == 2 || ord == 4) ;
    InMethod("RungeKutta::RungeKutta(ord)") ;
# endif

    order = ord ;
    k1    = NULL ;
    k2    = NULL ;
    k3    = NULL ;
    k4    = NULL ;
    uwork = NULL ;

    pastSolutions = new Vector<FlatField*>(1) ;
    pastSolutions->At(1) = NULL ;
}


RungeKutta :: ~RungeKutta ()
{
    // Destructor.

    unref(k1) ;
    unref(k2) ;
    unref(k3) ;
    unref(k4) ;
    unref(uwork) ;
}

FlatField* RungeKutta :: ComputeStrongHistory ()
{
    // Returns the "strong part" of the history vector at the current time step
    // n+1:
    //   Hist_strong(n+1) = u(n) / dt.
    // In order to save space, Hist_strong(n+1) overwrites u(n). Therefore this
    // method must be performed AFTER method ComputeWeakHistory, which uses u(n).

    pastSolutions->At(1)->Multiply(ONE/dt) ;

    return pastSolutions->At(1) ;
}


FlatField* RungeKutta :: ComputeWeakHistory ()
{
    // Returns the "weak part" of the history vector at the current time step
    // n+1: 
    //   if order = 2, returns: k2
    //   if order = 4, returns: (k1 + 2 k2 + 2 k3 + k4) / 6.
    // In order to save space, the answer overwrites k2 (if order=2) or k1 (if
    // order=4).
    // Programming note on how to compute, say, k3 from k2. The only trouble is
    //    to compute x = B-1 k2. This is done by solving the problem B x = k2,
    //    i.e., by solving the problem with lambda=1, mode=explicit, weak history
    //    =k2 and strong history=0.
    //    (in the case of a steady Stokes problem, this implies computing the
    //    pressure, e.g., inverting the Uzawa operator, and yields divergence-
    //    free x (and thus velocity) at every of the 2 or 4 sub-iterations)

    FlatField *uN, *answer ;
    real      tN ;

    uN = pastSolutions->At(1) ;
    tN = time - dt ;

    // settings for computing B-1 k1, B-1 k2 and B-1 k3
    semiDiscreteProblem->SetLambda(ONE) ; 
    semiDiscreteProblem->SetStrongHistory(NULL) ; 

    // k1
    semiDiscreteProblem->GetFunctional(uN,tN,k1) ;

    // k2
    semiDiscreteProblem->SetWeakHistory(k1) ; 
    semiDiscreteProblem->Solve() ;                  // B-1 k1 computed
    uwork->CopyFrom(uN) ;
    uwork->Add(solution->GetMain(),HALF*dt) ;       // uwork = uN + B-1 k1 dt/2
    semiDiscreteProblem->GetFunctional(uwork,tN+HALF*dt,k2) ;

    if (order == 2)
	answer = k2 ;

    else {  // order = 4 

	// k3 
	semiDiscreteProblem->SetWeakHistory(k2) ; 
	semiDiscreteProblem->Solve() ;                // B-1 k2 computed
	uwork->CopyFrom(uN) ;
	uwork->Add(solution->GetMain(),HALF*dt) ;     // uwork = uN + B-1 k2 dt/2
	semiDiscreteProblem->GetFunctional(uwork,tN+HALF*dt,k3) ;

	// k4
	semiDiscreteProblem->SetWeakHistory(k3) ; 
	semiDiscreteProblem->Solve() ;                // B-1 k3 computed
	uwork->CopyFrom(uN) ;
	uwork->Add(solution->GetMain(),dt) ;          // uwork = uN + B-1 k3 dt
	semiDiscreteProblem->GetFunctional(uwork,tN+dt,k4) ;

	// (k1+2k2+2k3+k4)/6
	answer = k1 ;
	answer->Add(k2,TWO) ;
	answer->Add(k3,TWO) ;
	answer->Add(k4) ;
	answer->Multiply(SIXTH) ;
    }

    semiDiscreteProblem->SetWeakHistory(NULL) ;
    return answer ;
}
    

void RungeKutta :: DoOneStep ()
{
    // Updates the temporal unknown by one time step, i.e., performs
    //   if order=2: u(n+1)/dt = u(n)/dt + k2 
    //   if order=4: u(n+1)/dt = u(n)/dt + (k1+2k2+2k3+k4)/6.

# ifdef REQUIRE
    Require("has a semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has initial values", solution != NULL) ;
    Require("has dt", dt != UNINITIALIZED_REAL) ;
    InMethod("RungeKutta::DoOneStep()") ;
# endif

    FlatField *strong, *weak ;
    printf(" RUNGE KUTTA DO ONE step \n");
 
    // housecleaning

    UpdateTime() ;

    // compute the history field Hist(n+1)
    weak   = ComputeWeakHistory() ;                        // overwrites k2 or k1
    strong = ComputeStrongHistory() ;                      // overwrites u(n)
    semiDiscreteProblem->SetWeakHistory(weak) ;
    semiDiscreteProblem->SetStrongHistory(strong) ;

    // set lambda (was overwritten in the computation of 'weak')
    semiDiscreteProblem->SetLambda(ONE/dt) ;

    // get u(n+1) ;
    semiDiscreteProblem->Solve() ;
}


void RungeKutta :: SetField (Field* x)
{
    // Sets to 'x' the starting values of the unknown field u.

    TimeIntegrationScheme::SetField(x) ;

    // if they exist, delete the k_i's, because their interpolation may differ
    // from that of 'x'.

    unref(uwork) ;
    unref(k1) ;
    unref(k2) ;
    unref(k3) ;
    unref(k4) ;

    // (re-)allocate the k_i's

    uwork = x->GetMain()->DuplicateEmpty() ;
    k1    = x->GetMain()->DuplicateEmpty() ;
    k2    = x->GetMain()->DuplicateEmpty() ;

    if (order == 4) {
	k3 = x->GetMain()->DuplicateEmpty() ;
	k4 = x->GetMain()->DuplicateEmpty() ;
    }
}


void RungeKutta :: SetSemiDiscreteProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the semi-discrete problem that the receiver will solve at
    // each time step.

# ifdef REQUIRE
    Require("has no semi-discrete problem yet", semiDiscreteProblem == NULL) ;
    Require("same mesh", solution == NULL || 
	    solution->GetMesh() == prob->GetMesh()) ;
    InMethod("RungeKutta::SetSemiDiscreteProblem(prob)") ;
# endif

    TimeIntegrationScheme::SetSemiDiscreteProblem(prob) ;

    semiDiscreteProblem->SetToExplicit() ;
}


void RungeKutta :: UpdateTime ()
{
    // Advances the date by one time increment and updates the past values.

    currentTimeStep++ ;
    time = time + dt ;
    semiDiscreteProblem->SetTime(time) ;

    // update 'pastSolution'
    pastSolutions->At(1)->CopyFrom(solution->GetMain()) ;
}
 


//_______________________________ ThetaMethod _________________________________


ThetaMethod :: ThetaMethod (real thet)
    : TimeIntegrationScheme ()
{
    // Constructor. Initializes the receiver to the rule with theta equal to 
    // 'thet'.

    theta = thet ;
    F     = NULL ;

    pastSolutions = new Vector<FlatField*>(1) ;
    pastSolutions->At(1) = NULL ;
}


ThetaMethod :: ~ThetaMethod ()
{
    // Destructor.

    unref(F) ;
}


FlatField* ThetaMethod :: ComputeStrongHistory ()
{
    // Returns the "strong" part of the history vector at the current time step
    // n+1:
    // - if theta != 0 (implicit): returns 1/(theta dt) u(n);
    // - if theta  = 0 (explicit): returns 1/dt u(n).
    // In order to save space, the answer overwrites u(n).

    FlatField *answer ;
    real      factor ;

    if (theta != ZERO)
	factor = ONE/(theta*dt) ;
    else
	factor = ONE/dt ;
    answer = pastSolutions->At(1) ;
    answer->Multiply(factor) ;

    return answer ;
}

  
FlatField* ThetaMethod :: ComputeWeakHistory ()
{
    // Returns the "weak" part of the history vector at the current time step
    // n+1:
    // - if theta != 0 (implicit): returns (1-theta)/theta F(n);
    // - if theta  = 0 (explicit): returns F(n).
    // In order to save space, the answer overwrites F(n).

    if (theta == ONE)
	return NULL ;

    if (theta != ZERO && theta != HALF)
	F->Multiply((ONE-theta)/theta) ;

    return F ;
}

  
void ThetaMethod :: DoOneStep ()
{
    // Performs one time step, i.e., solves for u(n+1):
    //   1/(theta dt) u(n+1) = F(u(n+1)) + Hist(n+1), if theta != 0 (implicit)
    //   1/dt u(n+1)         = Hist(n+1)              if theta  = 0 (explicit).
    // If the starting step has not been performed yet, performs it.

# ifdef REQUIRE
    Require("has a semi-discrete problem", semiDiscreteProblem != NULL) ;
    Require("has initial values", solution != NULL) ;
    Require("has dt", dt != UNINITIALIZED_REAL) ;
    InMethod("ThetaMethod::DoOneStep()") ;
# endif

    FlatField *weak, *strong ;

    // perform the starting step n=0 if necessary
    if (currentTimeStep == 0 && theta != ONE)
	semiDiscreteProblem->GetFunctional(solution->GetMain(),time,F) ;

    // housecleaning
    UpdateTime() ;

    // set lambda
    if (theta != ZERO)
	semiDiscreteProblem->SetLambda(ONE/(theta*dt)) ;
    else
	semiDiscreteProblem->SetLambda(ONE/dt) ;

    // compute the history fields
    weak   = ComputeWeakHistory() ;                 // overwrites F(n)
    strong = ComputeStrongHistory() ;               // overwrites u(n)
    semiDiscreteProblem->SetWeakHistory(weak) ;
    semiDiscreteProblem->SetStrongHistory(strong) ;

    // get u(n+1) ;
    semiDiscreteProblem->Solve() ;

    // compute F(u(n+1))
    if (theta != ONE)
	semiDiscreteProblem->GetFunctional(solution->GetMain(),time,F) ;
}


void ThetaMethod :: SetField (Field* x)
{
    // See description in version of base class TimeIntegrationScheme.

    TimeIntegrationScheme::SetField(x) ;

    // if they exist, delete 'F', because its interpolation may differ from that
    // of 'x'.
    unref(F) ;

    // (re-)allocate 'pastSolution' and 'F'
    if (theta != ONE)                          // F(n) not used in backward Euler
	F = x->GetMain()->DuplicateEmpty() ;
}


void ThetaMethod :: SetSemiDiscreteProblem (SemiDiscreteProblem* prob)
{
    // Sets to 'prob' the semi-discrete problem that the receiver will solve at
    // each time step.

# ifdef REQUIRE
    Require("has no semi-discrete problem yet", semiDiscreteProblem == NULL) ;
    Require("same mesh", solution == NULL || 
	    solution->GetMesh() == prob->GetMesh()) ;
    InMethod("ThetaMethod::SetSemiDiscreteProblem(prob)") ;
# endif

    TimeIntegrationScheme::SetSemiDiscreteProblem(prob) ;

    if (theta == ZERO)
	semiDiscreteProblem->SetToExplicit() ;
    else 
	semiDiscreteProblem->SetToImplicit() ;
}


void ThetaMethod :: UpdateTime ()
{
    // Advances the date by one time increment and updates the past values.

    currentTimeStep++ ;
    time = time + dt ;
    semiDiscreteProblem->SetTime(time) ;

    // update 'pastSolution'
    pastSolutions->At(1)->CopyFrom(solution->GetMain()) ;
}
 

