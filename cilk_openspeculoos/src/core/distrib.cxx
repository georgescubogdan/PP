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

// file distrib.cxx

#include "core/distrib.hxx"
#include "core/options.hxx"
#include <math.h>


//_________________________________ Distribution ______________________________


Distribution :: Distribution (int nelem)
    : RealVector (nelem+1), GarbagedObject()
{
    // Constructor. Initializes the receiver to a distribution with 'nelem'
    // segments.
}


int Distribution :: GetNbSegments () 
{
    // Returns the number of segments of the receiver.

    return size-1 ;
}


//_____________________________ UniformDistribution ___________________________


UniformDistribution :: UniformDistribution (int nelem)
    : Distribution (nelem)
{
    // Constructor.

# ifdef REQUIRE
    Require("'nelem' >= 1", nelem >= 1) ;
    InMethod("UniformDistribution::UniformDistribution(nelem)") ;
# endif

    real step ;
    int  i ;

    step = 2. / (GetSize()-1) ;
    for (i=0 ; i<GetSize()-1 ; i++) 
	values[i] = -1. + i*step ;
    values[GetSize()-1] = 1. ;              // avoid any round-off error
}


//___________________________ GeometricalDistribution _________________________


GeometricalDistribution :: GeometricalDistribution (int nelem, real ratio)
    : Distribution(nelem)
{
    // Constructor. Initializes the receiver to a geometrical distribution with
    // 'nelem' segments.
    // ratio = smallest size over largest size of intervals
    // the smallest interval is located near  1. if ratio is > 1
    //                                  near -1. if ratio is < 1
    // a and r are issued from:
    //     a + a r + a r^2 + a r^3 + ... + a r^(nelem-1) = 2
    //     a r^(nelem-1) / a = ratio 

    real step ;
    int  i ;

# ifdef REQUIRE
    Require("'nelem' >= 1", nelem >= 1) ;
    Require("'ratio' > 0.", ratio > 0.) ;
    InMethod("GeometricalDistribution::GeometricalDistribution(nelem)") ;
# endif

    if (nelem == 1) {
	values[0] = -1.; 
	values[1] =  1.;
    }
    else {
	r = exp(log(ratio) / (nelem-1)) ;
	if (fabs(r - 1.) > 1e-6)
	    a = 2.*(1.-r) / (1.-ratio*r) ; // normal case
	else 
	    a = 2. / nelem ;               // to avoid trouble if ratio is close to 1
	step      = a ;
	values[0] = -1. ;
	for (i=1 ; i<GetSize() ; i++) { 
	    values[i] = values[i-1] + step ;
	    step      = step * r ;
	}
	values[GetSize()-1] = 1. ;
    }
}

