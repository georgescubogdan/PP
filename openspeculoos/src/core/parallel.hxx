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

#ifndef PARALLEL_HEADER
#define PARALLEL_HEADER

/*!
 * \file parallel.hxx
 * \brief This file contains functions related to parallelism.
 * \author {Lots of authors}
 * \version 1.0
 *
 * This file contains functions related to parallelism.
 * For each function, two versions exist: one for a parallel implementation,
 * one for a scalar implementation. Whether the parallel or the scalar
 * version is compiled depends on the PARALLEL option. If the PARALLEL option
 * is used, you must compile and link with the communication librairies (e.g.,
 * mpich).
 *
 * Note: If PARALLEL is not defined, function Mpi_Wtime always returns 0.
 *       If PARALLEL is defined, function Mpi_Wtime seems to give the same
 *       result as an object of type Timer.
 *       In conclusion, using class Timer is simpler.
 *
 * Note: this file and its corresponding source file are heavily polluted by
 *       stuff related to "FCI" (i.e., to the Swiss-T0 project). If you are not
 *       interested in that project, you can remove everything related to it,
 *       thereby gaining some code readability.
 */

#include "core/vector.hxx"
#include "core/vector_templates.hxx"
#include "core/util.hxx"
#include "core/options.hxx"
#include <string>

#ifdef PARALLEL
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#endif


void    Mpi_Barrier () ;
void    Mpi_Communicate () ;
void    Mpi_Finalize () ;
void    Mpi_Initialize (int* argc, char*** argv) ;
inline boolean Mpi_IsOn (int p) ;
real    Mpi_Max (real x) ;
void    Mpi_PrintBuffers (int p=-1) ;
void    Mpi_Receive (int size, real* values, int incr, int proc) ;
void    Mpi_ReceiveInt (int* value, int proc) ;
void    Mpi_ReceiveVectInt (int size,int* value, int incr, int proc) ;
inline void    Mpi_ReceiveFromBuffer (int size, real* values, int incr, 
                                      int proc, int mode=0) ;
void    Mpi_ResetCommunication () ;
void    Mpi_Send (int size, real* values, int incr, int proc) ;
void    Mpi_SendInt (int* value, int proc) ;
inline void    Mpi_SendToBuffer (int size, real* values, int incr, int proc) ;
int     Mpi_Sum (int n) ;
real    Mpi_Sum (real x) ;
void Get_VecSum_all_proc(int* tabOut, int* tabIn, const int size) ;
void Get_VecSum_all_proc(real* tabOut, real* tabIn, const int size) ;
void Get_VecSum_all_proc(RealVector* tabOut, RealVector* tabIn, const int size) ;
real    Mpi_Wtime () ;

/*
 * Mpi_rank contains the rank (processor number) of the current processor.
 * On a scalar implementation its value is 0.
 */
extern int Mpi_rank ;


/*
 * Mpi_nbProcessors contains the number of processors.
 * On a scalar implementation its value is 1.
 */

extern int Mpi_nbProcessors ;


/*
 * . Mpi_sendBuffers[p] is the buffer that stores the values to be sent by the
 *   current processor to processor p during the next communication phase
 *   (in function Mpi_Communicate).
 * . Mpi_receiveBuffers[p] is the buffer that stores the values to be received 
 *   by the current processor from processor p during the next communication
 *   phase.
 * . Mpi_receiveBufferSizes[p] is the number of values to be received by the
 *   current processor from processor p during the next communication phase, 
 *   i.e., it is the size of Mpi_receiveBuffers[p].
 * On a scalar implementation, these lists are not used.
 */

extern Mpi_Buffer** Mpi_sendBuffers ;
extern Mpi_Buffer** Mpi_receiveBuffers ;
extern int*         Mpi_receiveBufferSizes ;


/*
 *__________________________________ inline __________________________________
 */


#ifdef PARALLEL


boolean Mpi_IsOn (int p)
{
/*
     * Returns 'true' if this method is invoked by processor 'p' (i.e., if the
     * rank of the current processor is 'p'), else returns 'false'.
 */
  
    return (Mpi_rank == p) ;
}


void Mpi_ReceiveFromBuffer (int size, real* values, int incr, int proc,
                            int mode)
{
/*
     * Receives 'size' coefficients from to buffer of processor 'proc', then:
     * - if 'mode'=0, copies these coefficients into 'values'
     * - if 'mode'=1, adds these coefficients to 'values'
     *
     * The initialized coefficients of 'values' are:
     *   values[0], values[incr], values[2*incr], etc.
     * The stride 'incr' can be negative.
 */

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid argument 'incr'", incr != 0) ;
    InMethod("Mpi_ReceiveFromBuffer(size,values,incr,proc)") ;
# endif

    (Mpi_receiveBuffers[proc])->Get(size,values,incr,mode) ;
}


void Mpi_SendToBuffer (int size, real* values, int incr, int proc)
{
/*
     * Sends to the buffer of processor 'proc' 'size' coeffs from 'values'.
     * The involved coeffs are:
     *   values[0], values[incr], values[2*incr], etc.
     * The stride 'incr' can be negative.
 */

# ifdef REQUIRE
    Require("involves 2 processors", Mpi_rank != proc) ;
    Require("valid argument 'incr'", incr != 0) ;
    InMethod("Mpi_SendToBuffer(size,values,incr,proc)") ;
# endif

    (Mpi_sendBuffers[proc])->Append(size,values,incr) ;
}


#else                                 // versions for a scalar implementation


boolean Mpi_IsOn (int p)
{
/*
     * Default implementation for a scalar computer.
     * Returns 'true' if 'p' = 0 (the number of the root processor), else returns
     * 'false'. Normally, it should always return 'true'.
 */
  
    return (p == 0) ;
}
  

void Mpi_ReceiveFromBuffer (int, real*, int, int, int)
{
/*
     * Default implementation for a scalar computer.
     * Issues an error message.
 */

    Error("Mpi_ReceiveFromBuffer(size,values,incr,proc,mode)",
	  "should not be called");
}


void Mpi_SendToBuffer (int, real*, int, int)
{
/*
     * Default implementation for a scalar computer.
     * Issues an error message.
 */

    Error("Mpi_SendToBuffer(size,values,incr,proc)","should not be called") ;
}


#endif       // end PARALLEL

#endif

