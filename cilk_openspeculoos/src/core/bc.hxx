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

#ifndef BOUNDARY
#define BOUNDARY

/*!
 * \file bc.hxx
 * \brief Boundary Conditions definitions
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/mesh.hxx"
#include "core/field.hxx"
#include "core/funct.hxx"
#include "core/vector.hxx"
#include "core/options.hxx"
#include <stdio.h>

#include "core/garbaged.hxx"
#include "core/param.hxx"

class Field ;
class FlatField ;
class Function ;
class Mesh ;
class RealVector ;


/*! \class DirichletCondition
   * \brief Manages Dirichlet boundary conditions.
   *
   * A Dirichlet condition is defined by:
   * - a prescribed-value field (attribute <{prescribedField}>),
   * - the number <{icomp}> of the component of the field it will apply to.
   * 
   * The prescribed values may be time-independent or time-dependent:
   * - in the former case, the prescribed values are given by attribute
   *  <{prescribedField}>.
   * - in the latter case, the user assigns a <{function}> f(x,y,z,t) to 
   * the Dirichlet condition. At any time t, the prescribed values are
   * given by attribute <{prescribedField}> with values initialized to
   * f(x,y,z,t).
   * 
   * An attribute <{homogeneousField}> is also needed. It is a copy of
   * <{prescribedField}>, but with values always initialized to 0. It is
   * used for applying effectively the condition in a homogeneous way.
   * 
   * A Dirichlet condition is applied to the icomp-th component of a target
   * field 'x' by means of method ApplyOn.
   * 
   * Caveat. Let dc be a Dirichlet condition to be applied on a field x.
   * It is important that the spectral elements of the mesh of attribute
   * <{prescribedField}> of dc be the same elements (i.e., the same 
   * objects, not just elements that are geometrically coincident) as the
   * elements of the mesh (or the mesh's skeleton) of x.
   * For example, if x is a 1D field and dc applies on a vertex v,
   * v must be pointed to by BOTH the mesh of x and the mesh of 
   * <{prescribedField}>. 
   * In pratice, when you get the mesh on which your Dirichlet condition is
   * to be created, get that mesh (or get its elements) using methods such
   * as GetSkeleton, GetSkin, ExtractMeshBetween, GetElementCloserTo, 
   * because these methods return existing elements, they never create
   * copies of existing elements.
   * If this condition is not respected, applying dc to x has no effect, 
   * because dc and x have no elements in common.
   */

class DirichletCondition : public GarbagedObject
{

protected :

    int        icomp ;
    FlatField* prescribedField ;
    Function*  function ;
    FlatField* homogeneousField ;

public :

// = CONSTRUCTORS

    DirichletCondition (FlatField* p, int ic) ;
    ~DirichletCondition () ;

// = SETTINGS

    void                SetFunction (Function* f) ;

// = APPLY BC

    void                ApplyOn (Field* x, real t=ZERO) ;
    void                ApplyHomogeneousOn (Field* x) ;
    void                ApplyHomogeneousOn (Field* x, int icx) ;

// = DUPLICATION

    DirichletCondition* DuplicateHomogeneous () ;

// = DEBUGGING

    void                Print () ;
} ;


/*! \class NeumannCondition
   * \brief Manages Neumann boundary conditions.
   *
   * A Neumann condition is defined by a prescribed flux field <{g}>. The
   * Neumann condition applies on the mesh of <{g}>.
   * 
   * The weak formulation of the Neumann condition involves:
   * - an <{integrationRule}> for performing the numerical quadratures. By
   * default, this integration rule is identical to (i.e, defined on the
   * same grid as) <{g}>.
   * - an external normal vector at every dof. Since the orientation of
   * that vector (external, as opposed to internal) cannot be defined
   * uniquely, an extra piece of information is stored: <{external}>. In
   * a 2D (resp. 3D) problem, <{external}> is a 2-component (resp.
   * 3-component) vector that gives (even approximately) the external
   * direction. At any dof, the sign of the normal is given by the 
   * following rule: the dot product of the normal and <{external}> must
   * be positive.
   * 
   * The user may want:
   * - to prescribe once and for all the values of <{g}>. In this case, 
   *   he/she just has to initialize <{g}>;
   * - to let the values of <{g}> vary, typically, to be time dependent. In
   *  this case, he/she supplies the Neumann condition with a
   * <{function}>. Whenever <{g}> is needed, the Neumann condition will
   *  initialize its values to <{function}>.
   * 
   * Finally, the problem on which the Neumann condition applies can be
   * defined on:
   * - a scalar (i.e., one-component) field, typically a pressure field,
   * - a vector (i.e., multi-component) field, typically a velocity field
   *   v.
   * In the latter case, the component of v on which <{g}> applies must be
   * specified: <{icomp}>.
   * 
   * Note that this scalar-vs-vector nature of the problem's field is
   * unrelated to the number of components of <{g}>, which must be equal to
   * the problem's dimension, e.g., 2 components in a 2D problem.
   */
class NeumannCondition : public GarbagedObject
{

protected :

    FlatField*  g ;
    int         icomp ;
    Function*   function ;
    RealVector* external ;
    FlatField*  integrationRule ;


public :

// = CONSTRUCTORS

    NeumannCondition (FlatField* grad, int ic, real ext1, real ext2, 
		      real ext3=ZERO) ;
    ~NeumannCondition () ;

// = SETTINGS

    void                SetIntegrationRule (FlatField* rule) ;
    void                SetFunction (Function* f) ;

// = APPLY BC

    void                ApplyOn (FlatField* x, real t=ZERO) ;
    FlatField*          GetIntegrationRule () ;

// = DEBUGGING
    void                Print (real time) ;
} ;

#endif

