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

#ifndef FUNCT
#define FUNCT

/*!
 * \file funct.hxx
 * \brief Functions
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/garbaged.hxx"
#include "core/param.hxx"

/*! \class Function
   * \brief Base class of classes managing user-supplied functions.
   * 
   * User-supplied functions are used to initialize fields, for example for
   * specifying exact values to a boundary-condition field.
   */ 

class Function : public GarbagedObject
{

public :

    virtual real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO,
                           int ic=1) ;
} ;

/*! \class Linear
   * \brief Implements linear functions
   * 
   * Implements linear functions :
   * - f(x)     = a + b x
   * - f(x,y)   = a + b x + c y
   * - f(x,y,z) = a + b x + c y + d z
   */ 
class Linear : public Function
{

protected :

    real aa ;
    real bb ;
    real cc ;
    real dd ;

public :

    Linear (real a, real b, real c=ZERO, real d=ZERO) ;
    
    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class Constant
   * \brief Implements constant functions
   * 
   * Implements constant functions:
   * - f(x)     = a
   * - f(x,y)   = a
   * - f(x,y,z) = a
   */ 

class Constant : public Linear
{

public :

    Constant (real a) ;
} ;


/*! \class QuadraticX
   * \brief Implements quadratic functions of x
   * 
   * Implements quadratic functions of x: f(x) = a x*x + b x + c.
   * 1D only: no 2D version f(x,y) nor 3D version f(x,y,z).
   */ 

class QuadraticX : public Function
{

protected :

    real aa ;
    real bb ;
    real cc ;

public :
 
    QuadraticX (real a, real b=ZERO, real c=ZERO) ;

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class QuadraticXY
   * \brief Implements quadratic functions
   * 
   * Implements quadratic functions:
   * f(x,y) = a1 x*x + a2 y*y + a3 x*y + a4 x + a5 y + a6
   * no f(x) nor f(x,y,z).
   */ 

class QuadraticXY : public Function
{

protected :

    real aa1 ;
    real aa2 ;
    real aa3 ;
    real aa4 ;
    real aa5 ;
    real aa6 ;

public :
 
    QuadraticXY (real a1, real a2, real a3, real a4, real a5, real a6) ;

    real Evaluate (real x, real y, real z=ZERO, real t=ZERO, int ic=1) ;
} ;



/*! \class Sin
   * \brief Implements sinus functions
   * 
   * Implements f(x)     = a sin(omega1 x)
   * or f(x,y)   = a sin(omega1 x) sin(omega2 y)
   * or f(x,y,z) = a sin(omega1 x) sin(omega2 y) sin(omega3 z),
   * depending on the number of 'omega' parameters transmitted to the
   * constructor.
   */ 
class Sin : public Function
{

protected :
 
    real aa ;
    real omega_1 ;
    real omega_2 ;
    real omega_3 ;

public :

    Sin (real a, real omega1, real omega2=UNINITIALIZED_REAL,
	 real omega3=UNINITIALIZED_REAL) ;

    real Evaluate (real x, real y=PI/TWO, real z=PI/TWO, real t=ZERO, 
                   int ic=1) ;
} ;


/*! \class Cos
   * \brief Implements cosinus functions
   * 
   * Implements f(x)     = a cos(omega1 x)
   * or f(x,y)   = a cos(omega1 x) cos(omega2 y)
   * or f(x,y,z) = a cos(omega1 x) cos(omega2 y) cos(omega3 z),
   * depending on the number of 'omega' parameters transmitted to the
   * constructor.
   */ 

class Cos : public Function
{

protected :
 
    real aa ;
    real omega_1 ;
    real omega_2 ;
    real omega_3 ;

public :

    Cos (real a, real omega1, real omega2=UNINITIALIZED_REAL,
	 real omega3=UNINITIALIZED_REAL) ;

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class CosSin
   * \brief Implements a cos(x) sin(y) functions
   * 
   * Implements f(x,y) = a cos(omega1 x) sin(omega2 y).
   */ 

class CosSin : public Function
{

protected :

    real aa ;
    real omega_1 ;
    real omega_2 ;

public :

    CosSin (real a, real omega1=ONE, real omega2=ONE) ;

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class SinCos
   * \brief Implements a sin(x) cos(y) functions
   * 
   * Implements f(x,y) = a sin(omega1 x) cos(omega2 y).
   */ 

class SinCos : public Function
{

protected :

    real aa ;
    real omega_1 ;
    real omega_2 ;

public :

    SinCos (real a, real omega1=ONE, real omega2=ONE) ;

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;

class SinExp : public Function
{
  // = TITLE
  //                                     b t
  //     Implements: f(x,t) = sin (a x) e   .
  //     Applies only in 1D.

  protected :

    real aa ;
    real bb ;

  public :

    SinExp (real a, real b) ;

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class RegularizedProfile
   * \brief Implements a regularized function
   * 
   * Implements a regularized function:
   * 
   * . if option = 0:
   *                                 m  2m
   *            f(x) = f_m(x) = (x (1-x))  2    on [0,1];
   * 
   * Note: f(0)=0, f(0.5)=1, f(1)=0.
   * 
   * . if option = 1:
   * 
   *            f(x,y) = f_m(x) f_n(y) = (x (1-x)) (y (1-y))  2    on [0,1]x[0,1]
   * 
   * Note: f(0,y) = f(1,y) = f(x,0) = f(x,1) = 0, f(0.5,0.5) = 1.
   * add by EP
   * . if option = 2
   * 
   *            f(x) = (1-x**2)**2  on [-1:1]
   */ 

class RegularizedProfile : public Function
{

protected :

    int  opt ;
    int  mm ;
    int  nn ;
    real power ;     //  2**(m+n), to speed up calculations

public :

    RegularizedProfile (int option, int m=0, int n=0) ;
    
    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class ExactHeat2D
   * \brief Implements a function used in a 2D heat-equation problem.
   * 
   * Implements a function used in a 2D heat-equation problem for giving
   * the exact expression of the unknown field:
   * 
   *                      3       -y
   *         f(x,y,t) = (x + 3) (e  + 1) sin(pi/2 t).
   * 
   *      See file main2D_heat.cxx.
   */ 

class ExactHeat2D : public Function
{

public :
   
    ExactHeat2D () ;

    real Evaluate (real x, real y, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class ForcingHeat2D
   * \brief Implements a function used in a 2D heat-equation problem.
   * 
   * Implements a function used in a 2D heat-equation problem for giving
   * the expression of the forcing term:
   * 
   *                           3      -y
   *         f(x,y,t) = pi/2 (x +3) (e  +1) cos(pi/2 t)
   * 
   *                          -y        3     -y
   *                  - (6x (e  +1) + (x +1) e  ) sin(pi/2 t)
   * 
   *      See file main2D_heat.cxx.
   */ 

class ForcingHeat2D : public Function
{

public :
   
    ForcingHeat2D () ;

    real Evaluate (real x, real y, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class ExactStokes2D
   * \brief Implements a function used in a 2D Stokes problem
   * 
   * Implements a function used in a 2D Stokes problem for giving the
   * exact expression of the velocity or of the pressure:
   * 
   *      - if option=1:
   * 
   *           f(x,y,t) = [ u(x,y,t) = -cx sy st ]
   *                      [ v(x,y,t) =  sx cy st ] 
   * 
   *      - if option=2:
   * 
   *           f(x,y,t) = p(x,y,t) = -pi sx sy st
   * 
   *      - if option=3:// ajout EP
   * 
   *           f(x,y,t) = [dp(x,y,t)/dx = -pi^2/2 cx sy st]
   *                      [dp(x,y,t)/dy = -pi^2/2 sx cy st]
   * 
   *      where cx stands for cos(pi/2 x), sx for sin(pi/2 x), etc.
   * 
   *      The same exact solution is also used in a 2D Navier-Stokes test (see
   *      file main2D_navier.cxx).
   * 
   *      Also used for a degenerated 3D case (main3D_stokes.cxx).
   */ 

class ExactStokes2D : public Function
{

protected :

    int option ;

public :
   
    ExactStokes2D (int opt) ;
    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;

/*! \class ExactPoiseuille2D
   * \brief Implements a function used in a 2D Stokes problem
   * 
   * Implements a function used in a 2D Stokes problem for giving the
   * exact expression of the velocity or of the pressure:
   */ 

class ExactPoiseuille2D : public Function
{

protected : 
    real aa ;
 
public :
   
    ExactPoiseuille2D (real a) ;
    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


class ForcingNavierStokes2D : public Function
{
  // = TITLE
  //                                     
  //     Implements a function used in a 2D Navier-Stokes problem for giving
  //     the exact expression of the forcing term of the momentum equation.
  //
  //     See file main2D_navier.cxx
  //
  //     The corresponding function for the solution (velocity and pressure) is
  //     class ExactStokes2D, the same as for class ForcingStokes2D.

  public :

    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;


/*! \class ForcingStokes2D
   * \brief Implements a function used in a 2D Stokes problem
   * 
   * Implements a function used in a 2D Stokes problem for giving the
   * exact expression of the forcing term of the momentum equation.
   * 
   * See file main2D_stokes.cxx.
   * 
   * Also used for a generated 3D case - see file main3D_stokes.cxx
   */ 

class ForcingStokes2D : public Function
{

public :
   
    real Evaluate (real x, real y=ZERO, real z=ZERO, real t=ZERO, int ic=1) ;
} ;



#endif
