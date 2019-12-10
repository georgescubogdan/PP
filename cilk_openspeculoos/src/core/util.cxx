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

// util.cxx

#include "core/util.hxx"
#include "core/point.hxx"
#include "core/parallel.hxx"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int nbViolated = 0 ;

// Initializes the variable 'nbViolated', which is the number of items in 
// 'violated'.


string violated ;

// 'violated' stores the messages to send to the user for every violated
// requirement (usually, 1).

real angl2pi(real COS,real SIN)
{
    real A;

    printf(" in a2pi cos = %lf sin =%lf \n", COS, SIN);
    if( (COS>=ZERO) && (SIN>=ZERO) && (COS < ONE)) 
	A = acos(COS);

    if((COS>=ZERO) && (SIN<ZERO) && (COS<ONE)) 
	A = TWO * PI - acos(COS);

    if((COS<ZERO) && (SIN>=ZERO) && (COS> ONE_MINUS ))  
	A = acos(COS);
    if((COS<ZERO) && (SIN<ZERO) && (COS > ONE_MINUS))   
	A = PI + acos(-COS);

    if((COS>=ONE))  { A = ZERO;}
    if((SIN>=ONE))  { A = HALF * PI  ;}
    if((COS<=ONE_MINUS)) { A = PI ;}
    if((SIN<=ONE_MINUS)) { A = 1.5 * PI;}

    printf(" end angle2pi A = %lf \n", A);
    return A;
}

int CheckConvexity (Point *bl, Point *br, Point *tl, Point *tr) 
{
    // In 2D (z=0), returns True if the 4 corners define a convex set.
    // 'bl', 'br', 'tl' and 'tr' stand for "bottom-left", "bottom-right", "top-
    // left" and "top-right", respectively.
    //
    // This is performed by evaluating the 4 possible cross product between adja-
    // cent vectors defined by corners. The quadrangle is convex if all those
    // products are of the same sign, provided that that they are evaluated in a
    // consistent way. 

    int answer = 1 ;

    Point b = *bl - *br ;
    Point r = *br - *tr ;
    Point t = *tr - *tl ;
    Point l = *tl - *bl ;

    answer = answer && ( (b % r).GetZ() > DEFAULT_TOLERANCE ) ;
    answer = answer && ( (r % t).GetZ() > DEFAULT_TOLERANCE ) ;
    answer = answer && ( (t % l).GetZ() > DEFAULT_TOLERANCE ) ;
    answer = answer && ( (l % b).GetZ() > DEFAULT_TOLERANCE ) ;

    return answer ;
}


void Error (string method, string message)
{
    // Issues an error message and interrupts the program execution.

    printf("\n\n******************* SPECULOOS ERROR ********************\n\n");
    printf("Error in %s :\n      %s.\n\n",method.c_str(),message.c_str()) ;
    printf("********************************************************\n\n\n") ;

    Exit(0) ;
}


real EvalLegendre (int n, real x)
{
    // Returns Pn(x), where Pn is the Legendre polynomial of degree 'n'.

# ifdef REQUIRE
    Require("'n' >= 0", n >= 0) ;
    InMethod("EvalLegendre(n,x)") ;
# endif

    int  k, m ;
    real sum ;

    m   = n/2 ;
    sum = ZERO ;
    for (k=0 ; k<=m ; k++)
	sum += Pow(ONE_MINUS,k) * Factorial(n+n-k-k) * Pow(x,n-k-k) /
	    (Pow(TWO,n) * Factorial(k) * Factorial(n-k) * Factorial(n-k-k)) ;

    return sum ;
}


void ExhaustedMemory ()
{
    // Issues an error message and interrupts the program execution.
    // This function is called by function set_new_handler, when operator "new"
    // is unable to deliver space for a new object, due to exhaustion of the free
    // store.
    // See Stroustrup, The C++ Programming Language, 2nd edition, p 99.

    printf("\n\n******************* SPECULOOS ERROR ********************\n\n");
    printf("Free store exhausted\n\n") ;
    printf("********************************************************\n\n\n") ;

    Exit(0) ;
}


void Exit (int code)
{
    // Calls the C routine exit().
    // Avoids including explicitly "stdlib.h", avoids compilation messages, 
    // avoids error messages in parallel execution mode.

    Mpi_Finalize() ;
    exit(code) ;
}


int Factorial (int n)
{
    // Returns n!.

# ifdef REQUIRE
    Require("n >= 0", n >= 0) ;
    InMethod("Factorial(n)") ;
# endif

    if (n <= 1)
	return 1 ;
    else
	return n * Factorial(n-1) ;
}


FILE* Fopen (const char* fileName, const char* mode)
{
    // Function 'fopen' if stdio.h, plus error checking.

    FILE *answer ;

    answer = fopen(fileName,mode) ;
    if (answer == NULL)
    {
	printf("File name = %s \n",fileName);
	Error("Fopen(fileName,mode)","cannot open the file") ;
    }
    return answer ;
}
 

void ImplementedBySubclasses (string method)
{
    // Issues an error message on the standard output saying that 'method' is a 
    // virtual method with no default implementation, then interrupts the program
    // execution.

    printf("\n\n******************* SPECULOOS ERROR ********************\n\n") ;
    if (Mpi_nbProcessors != 1)
	printf("                    on processor %d\n\n",Mpi_rank) ;
    printf("%s is implemented by subclasses.\n\n",method.c_str()) ;
    printf("********************************************************\n\n\n") ;
    Exit(0) ;
}


void InMethod (string method)
{
    // If the global variable 'violated' is empty, do nothing. If it contains
    // error messages, prints these on the standard output and interrupts the
    // program execution.

    int i ;

    if (nbViolated) {
	printf("\n\n******************* SPECULOOS ERROR ********************\n\n");
	if (Mpi_nbProcessors != 1)
	    printf("                    on processor %d\n\n",Mpi_rank) ;
	for (i=0 ; i<nbViolated ; i++)
	  printf("Violated requirement: %s\n",violated.c_str()) ;
	printf("In method:            %s\n\n",method.c_str()) ;
	printf("********************************************************\n\n\n") ;
	Exit(0) ;
    }
}


void NotImplemented (string method)
{
    // Issues an error message saying that 'method' is not implemented, then
    // interrupts the program execution.

    printf("\n\n******************* SPECULOOS ERROR ********************\n\n") ;
    if (Mpi_nbProcessors != 1)
	printf("                    on processor %d\n\n",Mpi_rank) ;
    printf("%s is not implemented.\n\n",method.c_str()) ;
    printf("********************************************************\n\n\n") ;

    Exit(0) ;
}


void NotImplemented (string method, int flag)
{
    // Issues an error message saying that some feature identified with a number
    // 'flag' is not implemented in 'method', then interrupts the program 
    // execution.

    printf("\n\n******************* SPECULOOS ERROR ********************\n\n") ;
    if (Mpi_nbProcessors != 1)
	printf("                    on processor %d\n\n",Mpi_rank) ;
    printf("In %s case %d is not implemented.\n\n",method.c_str(),flag) ;
    printf("********************************************************\n\n\n") ;

    Exit(0) ;
}


real ParentToPhysic (real x1, real x2, real r)
{
    // Mapping from r in [-1,1] to x in [x1,x2].

    return Phi1(r)*x1 + Phi2(r)*x2 ;
}


real Phi1 (real r)
{
    // Evaluates the linear function that is 1 at r=-1 and 0 at r=1.

    return (ONE-r)*HALF ;
}


real Phi2 (real r)
{
    // Evaluates the linear function that is 0 at r=-1 and 1 at r=1.

    return (ONE+r)*HALF ;
}


real Pow (real x, int n)
{
    // Returns x**n.
    // Returns 1 if x=0 and n=0.

    if (x==ZERO && n==0)
	return ONE ;
    else
	return pow(x,(real)n) ;
}


void Require (string message, boolean test)
{
    // If 'test' is equal to False, stores 'message' in the list of violated
    // requirements. This list will be printed and the program execution
    // interrupted next time function 'InMethod' is invoked.

    string s ;

    if (! test) {
      s = message ;
	violated = s ;
    }
} 


void Warning (string method, string message)
{
    // Issues a warning message to the user, on the standard output.
    // Execution is not interrupted.

    printf("\n\n****************** SPECULOOS WARNING *******************\n\n") ;
    if (Mpi_nbProcessors != 1)
	printf("                    on processor %d\n\n",Mpi_rank) ;
    printf("Warning in %s : \n  %s\n\n",method.c_str(),message.c_str()) ;
    printf("********************************************************\n\n\n") ;
}


