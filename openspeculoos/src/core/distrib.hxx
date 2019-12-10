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

#ifndef DISTRIB
#define DISTRIB
 

/*!
 * \file distrib.hxx
 * \brief Distribution of parameters upon interval [-1,1]
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/garbaged.hxx"
#include "core/param.hxx"
#include "core/vector.hxx"


/*! \class Distribution
   * \brief Manages distributions of a parameter on the interval [-1,1].
   * 
   * The number of segments (elements) is required for the construction of 
   * an object.
   */
class Distribution : public RealVector, public GarbagedObject
{

public :

    Distribution (int nelem) ;

    int GetNbSegments() ;
} ;


/*! \class UniformDistribution
   * \brief Manages uniform distributions of points on the interval [-1,1].
   */

class UniformDistribution : public Distribution 
{

public :

    UniformDistribution (int) ;
} ;


/*! \class GeometricalDistribution
   * \brief Manages geometrical distributions of points on the interval [-1,1].
   * 
   * A geometrical distribution is defined by means of a 'reason' <{a}> and
   * a 'ratio' <{r}> such that
   * a + a r + a r^2 + a r^3 + ... + a r^(nelem-1) = 2
   */
class GeometricalDistribution : public Distribution 
{

protected :

    real a ;
    real r ;

public :
   
    GeometricalDistribution (int nelem, real ratio) ;
} ;


#endif
