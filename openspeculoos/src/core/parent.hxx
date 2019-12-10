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

#ifndef PARENT
#define PARENT

/*!
 * \file parent.hxx
 * \brief Manages parent element types
 * \author {Lots of authors}
 * \version 1.0
 */

#include "core/vector.hxx"
#include "core/vector_templates.hxx"
#include <string>

class ListMatrices ;
class Matrix ;
class FullMatrix ;

enum ParentElementType { GLL_PT, GLLI_PT, IGLL_PT, GL_PT, HP_PT} ;


/*! \class ParentElement
   * \brief Base class managing parent element types.
   * 
   * A parent element is a mono-dimensional information that is intended to
   * provide tools for interpolations and interpolated derivatives on the
   * reference interval [-1,1].
   * When the interpolation type is associated to a quadrature rule, the
   * quadrature weights are also provided by this class.
   * A derived class defines a new parent element type. Inside this class,
   * parent-element objects differ by their polynomial degree.
   * Parent elements are combined (one per direction) to define 
   * interpolations on spectral elements.
   * 
   * Two parent elements of the same kind and of the same degree have 
   * identical collocation points, quadrature weights, interpolation and 
   * derivative matrices. Therefore, in order to avoid duplications, parent
   * elements are internally arranged in vectors, one (static) vector per
   * derived class (<{instances}>). The i-th component of this vector is 
   * the parent element of degree i. This parent element can be pointed to
   * by several elementary fields, or several times by the same elementary
   * field.
   * In order to get a parent element, say the GLL parent of degree i, the
   * message "ParentElementGLL::GetParentElementOfDegree(i)" should be 
   * used, rather than a call to the constructor ParentElementGLL(i).
   * 
   * Computing the interpolation matrix (or the interpolated-derivative
   * matrix) from a parent element to another one is time consuming.
   * Therefore every parent element stores lists (<{listInterpMatrices}>
   * and <{listInterpDerivMatrices}>) of interpolation matrices to the
   * other parent elements. Every such matrix H can be used for computing
   * products x_target = H * x_source, where x_source is a vector of values
   * on the owner of M.
   */ 

class ParentElement 
{

protected :

    int           polynomialDegree ;
    RealVector*   collocationPoints ;
    RealVector*   quadratureWeights ;
    FullMatrix*   derivMatrix ;
    ListMatrices* listInterpMatrices ;
    ListMatrices* listInterpDerivMatrices ;

public :

    // = CONSTRUCTORS

    ParentElement (int n) ;
    virtual ~ParentElement () ;

    // = MANAGEMENT OF INSTANCES

    ParentElement*            GetParentElementOfDegree_ (int n) ;
    static void               DeleteAllInstances () ;

    // = INQUIRY

    virtual string            GetName (string s) = 0 ;
    inline  int               GetPolynomialDegree () ;
    inline  int               GetNbCollocationPoints () ;
    virtual int               GetNbInteriorCollocationPoints () = 0 ;
    virtual ParentElementType GetType () = 0 ;
    virtual boolean           IsDefinedOnBorders () ;

    // = QUADRATURE

    virtual RealVector*       GetWeights () ;
    inline RealVector*        GetCollocationPoints () ;

    // = DERIVATION
   
    virtual FullMatrix*       GetDerivMatrix () = 0 ;

    // = INTERPOLATION, INTERPOLATION OF DERIVATIVE

    virtual Matrix*           GetInterpMatrixTo (ParentElement* p) = 0 ;
    virtual FullMatrix*       GetInterpDerivMatrixTo (ParentElement* p) ;

    // = DEBUGGING

    void                      Print () ; 
    FullMatrix*               ComputeInterpMatrixTo (ParentElement* p) ;
    FullMatrix*               ComputeInterpDerivMatrixTo (ParentElement* p) ;
    real                      EvaluatePolynomial (int j, real x) ;
    real                      EvaluatePolynomialDerivative (int j, real x) ;
} ;

/*! \class ParentElementGLL
   * \brief Manages Gauss-Lobatto-Legendre parent elements.
   * 
   * Collocation points, quadrature weights, interpolation and derivative
   * matrices are computed by means of the MIT FORTRAN library developed by
   * E. Ronquist.
   */ 
class ParentElementGLL : public ParentElement
{ 

protected :

    static Vector<ParentElementGLL*>* instances ;

    ParentElementGLL (int n) ;
    ~ParentElementGLL () ;

public :

    // = MANAGEMENT OF INSTANCES

    static ParentElement*    GetParentElementOfDegree (int n) ;
    static void              DeleteAllInstances () ;

    // = INQUIRY
    
    string                   GetName (string s) ;
    inline int               GetNbInteriorCollocationPoints () ;
    inline ParentElementType GetType () ;
    boolean                  IsDefinedOnBorders () ;         

    // = DERIVATION

    FullMatrix*              GetDerivMatrix () ;

    // = INTERPOLATION

    Matrix*                  GetInterpMatrixTo (ParentElement* p) ;
} ;



/*! \class ParentElementGLLI
   * \brief Manages interior Gauss-Lobatto-Legendre parent elements.
   * 
   * Removes the external points of the GLL parent elements.
   * Uses the same weights as the GLL parent elements for quadrature.
   */ 

class ParentElementGLLI : public ParentElement
{

protected :

    RealVector*   collocationPointsGLL ;
    static Vector<ParentElementGLLI*>* instances ;

    ParentElementGLLI (int n) ;
    ~ParentElementGLLI () ;

public :
   
    // = MANAGEMENT OF INSTANCES

    static ParentElement*    GetParentElementOfDegree (int n) ;
    static void              DeleteAllInstances () ;

    // = INQUIRY
    
    string                   GetName (string s) ;
    inline int               GetNbInteriorCollocationPoints () ;
    inline ParentElementType GetType () ;
//    RealVector*              GetWeights () ;

    // = DERIVATION

    FullMatrix*              GetDerivMatrix () ;

    // = INTERPOLATION

    Matrix*                  GetInterpMatrixTo (ParentElement* p) ;
} ;


/*! \class ParentElementIGLL
   * \brief Manages interior Gauss-Lobatto-Legendre parent elements.
   * 
   * Removes the external points of the GLL parent elements.
   * Interpolation is provided but not derivation.
   * 
   * FIXME : [Vincent Keller] What's de difference with the previous one ?
   */

class ParentElementIGLL : public ParentElement
{

protected :

    static Vector<ParentElementIGLL*>* instances ;

    ParentElementIGLL (int n) ;
    ~ParentElementIGLL () ;

public :
   
    // = MANAGEMENT OF INSTANCES

    static ParentElement*    GetParentElementOfDegree (int n) ;
    static void              DeleteAllInstances () ;

    // = INQUIRY

    string                   GetName (string s) ;
    inline int               GetNbInteriorCollocationPoints () ;
    inline ParentElementType GetType () ;
    RealVector*              GetWeights () ;

    // = DERIVATION

    FullMatrix* GetDerivMatrix () ;

    // = INTERPOLATION

    Matrix*                  GetInterpMatrixTo (ParentElement* p) ;
} ;


/*! \class ParentElementGL
   * \brief Manages Gauss-Legendre parent elements.
   * 
   * Collocation points, quadrature weights, interpolation and derivative 
   * matrices are computed by means of the MIT FORTRAN library developed by
   * E. Ronquist.
   */ 
class ParentElementGL : public ParentElement
{

protected :

    static Vector<ParentElementGL*>* instances ;

    ParentElementGL (int n) ;
    ~ParentElementGL () ;

public :

    // = MANAGEMENT OF INSTANCES

    static ParentElement*    GetParentElementOfDegree (int n) ;
    static void              DeleteAllInstances () ;

    // = INQUIRY
    
    string                   GetName (string s) ;
    inline int               GetNbInteriorCollocationPoints () ;
    inline ParentElementType GetType () ;

    // = DERIVATIVE

    FullMatrix*              GetDerivMatrix () ;

    // = INTERPOLATION

    Matrix*                  GetInterpMatrixTo (ParentElement* p) ;
} ;


/*! \class ParentElementHP
   * \brief Manages h-p parent elements.
   * 
   * This class implements only the interpolation matrix and the
   * interpolated derivative matrix.
   * Collocation points, quadrature weights and derivative matrix are not
   * relevant.
   */ 
class ParentElementHP : public ParentElement
{ 

protected :

    static Vector<ParentElementHP*>* instances ;

    ParentElementHP (int) ;                    // protected methods
    ~ParentElementHP () ;
    real EvalHP (int j, real r) ;
    real EvalDerivHP (int j, real r) ;

public :

    // = MANAGEMENT OF INSTANCES

    static ParentElement*    GetParentElementOfDegree (int p) ;
    static void              DeleteAllInstances () ;

    // = INQUIRY

    string                   GetName (string s) ;
    inline int               GetNbInteriorCollocationPoints () ;
    inline ParentElementType GetType () ;
    boolean                  IsDefinedOnBorders () ;         

    // = DERIVATION

    FullMatrix*              GetInterpDerivMatrixTo (ParentElement* p) ;
    FullMatrix*              GetDerivMatrix () ;

    // = INTERPOLATION

    Matrix*                  GetInterpMatrixTo (ParentElement*) ;
} ;


//_____________________________ ParentElement (inline) ________________________


RealVector* ParentElement :: GetCollocationPoints() 
{
    // Returns the array of the positions (between -1 and 1) of the collocation
    // points of the receiver.

    return collocationPoints ; 
}


int ParentElement :: GetNbCollocationPoints ()
{
    // Returns the number of collocation points of the receiver.
    // Programming note: the choice of putting a test on the nature of the
    //   receiver (h-p or not h-p) instead of making the function virtual was
    //   dictated by efficiency considerations.

# ifdef REQUIRE
    Require("not an h-p parent", GetType() != HP_PT) ;
    InMethod("ParentElement::GetNbCollocationPoints ()") ;
# endif

    return polynomialDegree + 1 ;
}


int ParentElement :: GetPolynomialDegree() 
{
    // Returns the polynomial degree of the receiver.

    return polynomialDegree ; 
}


//___________________________ ParentElementGLL (inline) _______________________


int ParentElementGLL :: GetNbInteriorCollocationPoints() 
{ 
    return max(polynomialDegree-1,0) ; 
}


ParentElementType ParentElementGLL :: GetType() 
{ 
    return GLL_PT ; 
}


//___________________________ ParentElementGLLI (inline) ______________________


int ParentElementGLLI :: GetNbInteriorCollocationPoints() 
{ 
    return polynomialDegree + 1 ; 
}


ParentElementType ParentElementGLLI :: GetType() 
{ 
    return GLLI_PT ; 
}


//___________________________ ParentElementIGLL (inline) ______________________


int ParentElementIGLL :: GetNbInteriorCollocationPoints() 
{ 
    return polynomialDegree + 1 ; 
}


ParentElementType ParentElementIGLL :: GetType() 
{ 
    return IGLL_PT ; 
}


//__________________________ ParentElementGL (inline) _________________________


int ParentElementGL :: GetNbInteriorCollocationPoints () 
{ 
    return polynomialDegree + 1 ; 
}


ParentElementType ParentElementGL :: GetType () 
{ 
    return GL_PT ;
}


//__________________________ ParentElementHP (inline) _________________________


int ParentElementHP :: GetNbInteriorCollocationPoints () 
{
    // Issues an error message, since h-p basis are not localized.

    NotImplemented("ParentElementHP::GetNbInteriorCollocationPoints()") ;
    return 0 ;
}


ParentElementType ParentElementHP :: GetType () 
{ 
    return HP_PT ;
}


#endif

