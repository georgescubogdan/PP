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

#ifndef MAINHEADERS
#define MAINHEADERS

/*!
 * \file main.hxx
 * \brief Main includes file
 * \author {Lots of authors}
 * \version 1.0
 */

/*! \mainpage OpenSpecuLOOS 
 *
 * FIXME: This is just a bulk !!!
 *
 * \section intro_sec Introduction
 *
 * FIXME: This is the introduction.
 *
 * \section install_sec Installation
 *
 * For more information on the Installation, see INSTALL file at the root of the distribution.
 *
 * \subsection step1 Step 1: Create a main.cxx file.
 *
 * Depending on which problem you want to solve.
 *
 * \subsection step2 Step 2: Create a link
 *
 * In the /main directory, create a symbolic link using 
 * ln -s main.cxx your_main_file.cxx
 * 
 * \subsection step3 Step 1: Compile
 *  
 * Sequence ./configure, make, make install
 *  
 * \section help Help
 * 
 * Contact us for more information.
 * 
 */


#include "core/bc.hxx"
#include "core/communic.hxx"
#include "core/distrib.hxx"
#include "core/edge.hxx"
#include "core/element.hxx"
#include "core/face.hxx"
#include "core/field.hxx"
#include "core/funct.hxx"
#include "garbaged.hxx"
#include "jacobi.hxx"
#include "matrix.hxx"
#include "mesh.hxx"
#include "parallel.hxx"
#include "parent.hxx"
#include "point.hxx"
#include "postpro.hxx"
#include "precond.hxx"
#include "problem.hxx"
#include "solver.hxx"
#include "timer.hxx"
#include "timinteg.hxx"
#include "vertex.hxx"
#include "vector.hxx"
#include "volume.hxx"
#include "timer.hxx"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include <unistd.h>



using namespace std ;

#endif
