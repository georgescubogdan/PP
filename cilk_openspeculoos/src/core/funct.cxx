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

// function.cxx

#include "core/funct.hxx"
#include "core/util.hxx"
#include <math.h>


//__________________________________ Function _________________________________


real Function :: Evaluate (real, real, real, real, int)
{
    // Issues an error message.

    NotImplemented("Function::Evaluate(x,y,z,t)") ;
    return 0. ;
}


//___________________________________ Linear __________________________________


Linear :: Linear (real a, real b, real c, real d)
{
    // Constructor.

    aa = a ;
    bb = b ;
    cc = c ;
    dd = d ;
}


real Linear :: Evaluate (real x, real y, real z, real, int)
{
    return aa + bb*x + cc*y + dd*z ;
}


//__________________________________ Constant _________________________________


Constant :: Constant (real a)
    : Linear(a,ZERO,ZERO,ZERO) 
{
    // Constructor.
}


//________________________________ QuadraticX _________________________________


QuadraticX :: QuadraticX (real a, real b, real c)
{
    // Constructor.

    aa = a ;
    bb = b ;
    cc = c ;
}


real QuadraticX :: Evaluate (real x, real y, real z, real, int)
{
    Require("y = z = 0", y == 0. && z == 0.) ;
    InMethod("QuadraticX::Evaluate(x,y,z)") ;

    return aa*x*x + bb*x + cc ;
}


//________________________________ QuadraticXY ________________________________


QuadraticXY :: QuadraticXY (real a1, real a2, real a3, 
                            real a4, real a5, real a6)
{
    // Constructor.

    aa1 = a1 ;
    aa2 = a2 ;
    aa3 = a3 ;
    aa4 = a4 ;
    aa5 = a5 ;
    aa6 = a6 ;
}


real QuadraticXY :: Evaluate (real x, real y, real z, real, int)
{
    Require("z = 0", z == 0.) ;
    InMethod("QuadraticXY::Evaluate(x,y,z)") ;

    return aa1*x*x + aa2*y*y + aa3*x*y + aa4*x + aa5*y + aa6 ;
}


//____________________________________ Sin ____________________________________


Sin :: Sin (real a, real omega1, real omega2, real omega3)
{
    // Constructor.

    aa      = a ;
    omega_1 = omega1 ;
    omega_2 = omega2 ;
    omega_3 = omega3 ;
}


real Sin :: Evaluate (real x, real y, real z, real, int)
{
    real answer ;

    answer = aa * sin(omega_1*x) ;

    if (omega_2 != UNINITIALIZED_REAL) 
        answer *= sin(omega_2*y) ;

    if (omega_3 != UNINITIALIZED_REAL) 
        answer *= sin(omega_3*z) ;

    return answer ;
}


//____________________________________ Cos ____________________________________


Cos :: Cos (real a, real omega1, real omega2, real omega3)
{
    // Constructor.

    aa      = a ;
    omega_1 = omega1 ;
    omega_2 = omega2 ;
    omega_3 = omega3 ;
}


real Cos :: Evaluate (real x, real y, real z, real, int)
{
    real answer ;

    answer = aa * cos(omega_1*x) ;

    if (omega_2 != UNINITIALIZED_REAL) 
        answer *= cos(omega_2*y) ;

    if (omega_3 != UNINITIALIZED_REAL) 
        answer *= cos(omega_3*z) ;

    return answer ;
}


//__________________________________ CosSin ___________________________________


CosSin :: CosSin (real a, real omega1, real omega2)
{
    // Constructor.

    aa      = a ;
    omega_1 = omega1 ;
    omega_2 = omega2 ;
}


real CosSin :: Evaluate (real x, real y, real, real, int)
{
    return aa * cos(omega_1*x) * sin(omega_2*y) ;
}


//__________________________________ SinCos ___________________________________


SinCos :: SinCos (real a, real omega1, real omega2)
{
    // Constructor.

    aa      = a ;
    omega_1 = omega1 ;
    omega_2 = omega2 ;
}


real SinCos :: Evaluate (real x, real y, real, real, int)
{
    return aa * sin(omega_1*x) * cos(omega_2*y) ;
}

//__________________________________ SinExp ___________________________________


SinExp :: SinExp (real a, real b)
{
  aa = a ;
  bb = b ;
}


real SinExp :: Evaluate (real x, real, real, real t, int)
{
  return sin (aa*x) * exp(bb*t) ;
}

//_______________________________ RegularizedProfile __________________________


RegularizedProfile :: RegularizedProfile (int option, int m, int n)
{
    opt   = option ;
    mm    = m ;
    nn    = n ;
    power = pow(TWO,m+m+n+n) ;
}


real RegularizedProfile :: Evaluate (real x, real y, real z, real, int)
{
    if (opt == 0)
        return pow(x*(ONE-x),mm) * power ;
    else  if (opt == 1)
        return (ONE-pow(TWO*(x-HALF),mm)) * (ONE-pow(TWO*(x-HALF),mm)) * 
               (ONE-pow(TWO*(z-HALF),nn)) * (ONE-pow(TWO*(z-HALF),nn)) ;
    else if (opt == 2) // option == 2
        return ((ONE - x * x) * (ONE - x * x));
    else // if (opt == 3) // regularised profil on [0,1]
        return 16 * x * x * (ONE - x ) * (ONE - x );
}


//___________________________ Poiseuille 2D __________________________________

   
ExactPoiseuille2D :: ExactPoiseuille2D (real a) 
{
    aa = a ;
}

real ExactPoiseuille2D :: Evaluate (real x, real y, real, real, int) {
    return aa*(y*y -ONE)+0.*x ;
} 


//________________________________ ExactHeat2D ________________________________


ExactHeat2D :: ExactHeat2D ()
{
}


real ExactHeat2D :: Evaluate (real x, real y, real, real t, int)
{
    return (x*x*x + THREE) * (exp(-y) + ONE) * sin (HALF * PI * t) ;
}


//____________________________ ForcingNavierStokes2D __________________________


real ForcingNavierStokes2D :: Evaluate (real x, real y, real, real t, int ic)
{
  real cx, cy, ct, sx, sy, st, c1, c2, u, v, ux, uy, vx, vy ;
  real P2 = PI/TWO ;

  cx = cos(P2*x) ;
  cy = cos(P2*y) ;
  ct = cos(P2*t) ;
  sx = sin(P2*x) ;
  sy = sin(P2*y) ;
  st = sin(P2*t) ;

  u  = -cx*sy*st ;                              // exact velocity u
  v  =  sx*cy*st ;                              // exact velocity v

  if (ic == 1) {                                // x component 
    c1 = -PI*PI*cx*sy*st - P2*cx*sy*ct ;        //   Stokes contribution
    ux =  sx*sy*st*P2 ;
    uy = -cx*cy*st*P2 ;
    c2 = u*ux + v*uy ;                          //   convective contribution
  }
  else if (ic==2) {                                        // y component
    c1 = P2*sx*cy*ct ;                          //   Stokes contribution
    vx =  cx*cy*st*P2 ;
    vy = -sx*sy*st*P2 ;
    c2 = u*vx + v*vy ;                          //   convective contribution
  }
  else {
      c1 = ZERO ;
      c2 = ZERO ;
  }

  return c1+c2 ;
}


//_______________________________ ForcingHeat2D _______________________________


ForcingHeat2D :: ForcingHeat2D ()
{
}


real ForcingHeat2D :: Evaluate (real x, real y, real, real t, int)
{
    real ey, x3p3, pi2t ;

    ey   = exp(-y) ;
    x3p3 = x*x*x + THREE ;
    pi2t = HALF * PI * t ;

    return  HALF  * PI * x3p3 * (ey+ONE) * cos(pi2t)
        - (SIX*x*(ey+ONE) + x3p3*ey) * sin(pi2t) ;
}


//_______________________________ ExactStokes2D _______________________________


ExactStokes2D :: ExactStokes2D (int opt)
{
    option = opt ;
}


real ExactStokes2D :: Evaluate (real x, real y, real, real t, int ic)
{
    real P2 = PI/TWO ;

    if (option == 1) {        // velocity
        if (ic == 1)
            return -cos(P2*x) * sin(P2*y) * sin(P2*t) ;
        else if (ic == 2)
            return  sin(P2*x) * cos(P2*y) * sin(P2*t) ;
        else
            return 1.5 ;
    }
    else if (option == 2)                    // pressure
        return -PI * sin(P2*x) * sin(P2*y) * sin(P2*t) ;

    else if (option == 3){ //  (option == 3) {        // grad p
        if (ic == 1)
            return -PI*P2*cos(P2*x) * sin(P2*y) * sin(P2*t) ;
        else if (ic == 2)
            return -PI*P2*sin(P2*x) * cos(P2*y) * sin(P2*t) ;
        else
            return 1.5 ;
    }
    else { //  (option == 4) {        // dv/dt
        if (ic == 1)
            return -P2*cos(P2*x) * sin(P2*y) * cos(P2*t) ;
        else if (ic == 2)
            return  P2*sin(P2*x) * cos(P2*y) * cos(P2*t) ;
        else
            return 1.5 ;
    }

}


    

//______________________________ ForcingStokes2D ______________________________


real ForcingStokes2D :: Evaluate (real x, real y, real, real t, int ic)
{
    real cx, cy, ct, sx, sy, st ;
    real P2 = PI/TWO ;

    cx = cos(P2*x) ;
    cy = cos(P2*y) ;
    ct = cos(P2*t) ;
    sx = sin(P2*x) ;
    sy = sin(P2*y) ;
    st = sin(P2*t) ;

    if (ic == 1)
        return - PI*PI*cx*sy*st - P2*cx*sy*ct ;
    else if (ic == 2)
        return + P2*sx*cy*ct ;
    else
        return 0. ;
}
  
