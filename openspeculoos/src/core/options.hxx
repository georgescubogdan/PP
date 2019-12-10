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

#ifndef OPTIONS
#define OPTIONS

/*!
 *  \file options.hxx
 *  \brief This file defines several compilation options.
 *  \author {Lots of authors}
 *  \version 1.0
 * 1) The REQUIRE option:
 * ----------------------
 *  Defining the REQUIRE option implies that at the entry point of many methods
 * run-time checking will be performed on the arguments (and sometimes on the
 * object's attributes). Typically range-checking is performed when using arrays
 * or matrices.
 * In a development phase, it is advisable to set REQUIRE to: defined. Undefine
 * it when you really strive for (or measure) CPU-time efficiency.

 * 2) The CHECK option:
 * --------------------
 * Defining the CHECK option implies that with certain methods checks will be
 * performed on potentially wrong variables. 
 * In a development phase, it is advisable to set CHECK to: defined. Undefine
 * it when you really strive for (or measure) CPU-time efficiency.

 * 3) The PARALLEL option:
 * -----------------------
 * Defining the PARALLEL option generates a parallel executable, i.e., an
 * executable for a multiprocessor machine. Processor communication is ensured by
 * calls to the MPI (Message Passing Interface) library.
 * Not defining PARALLEL generates a scalar executable, i.e., an executable for
 * a scalar computer or for a 1-processor execution on a parallel machine.
 * If you define PARALLEL you must make sure that an implementation of the MPI
 * library is installed on your computer; you may also have to modify the paths
 * to the MPI library in the makefile.
 * By default, do not define PARALLEL.

 * 4) The FAST_GRADIENTS option:
 * -----------------------------
 * Defining the FAST_GRADIENTS options implies that the computations of 
 * gradients will be speeded up for straight elements, i.e., for elements whose
 * inverse jacobian matrix contains zeroes.
 * Undefine it typically if you suspect that the small perturbation this option 
 * implies is a cause of non-convergence in conjugate-gradient iterations.
 * By default, define FAST_GRADIENTS.

 */


//#define REQUIRE
//#define CHECK
#define FAST_GRADIENTS
#define FILTER

#endif
