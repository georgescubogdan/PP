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

#ifndef TIMER
#define TIMER

/*!
 * \file timer.hxx
 * \brief Clock Timer functionnalities
 * \author {Lots of authors}
 * \version 1.0
 */


#include "core/param.hxx"
#include <stdio.h>
#include <sys/times.h>
#include <fstream>
#include <string>
#include <map>
#include "core/parallel.hxx"

/*! \class Timer
   * \brief Provides clock timer functionalities (user, system, elapsed).
   * 
   * Measured timings is provided by calls to methods <{start}> and
   * <{stop}>. 
   * A timer can be started and stopped several times, yielding accumulated
   * times.
   * Relies on the <{times}> C function.
   */ 

class Timer
{

protected :

    string   nameTimer ;
    float    cpu_clktck ;
    boolean  isActive ;
    long     nb_calls ;
    tms      start_cstime ;
    tms      stop_cstime ;
    tms      tot_cstime ;
    clock_t  start_etime ;
    clock_t  stop_etime ;
    clock_t  tot_etime ;

public :

    Timer () ;
    Timer (string aname) ;
    ~Timer () ;

    void Start () ; 
    void Stop () ;
    void Print () ;
    void PrintIn (FILE* file) ;
    void GetTimeFloat (float&, float&,float& );
    float GetCPU ();

    // Add-on Vincent Keller
    void PrintSeconds () ;


} ;

/*! \class ParallelTimer
   * \brief Provides clock timer functionalities (user, system, elapsed) in parallel environment
   * 
   * The ParallelTimer is based on the getTimeOfDay() C method and should work only 
   * on UNIX based systems (FIXME : [Vincent Keller] verification ....
   */ 

#ifdef PARALLEL

class ParallelTimer{

    double start_time;
    double stop_time;
    double inter_time;
    
public :
    
    ParallelTimer () ;
    ~ParallelTimer () ;

    // timer related subs
	
    void Start () ; 
    void Stop () ;
    double GetIntermediateTime () ;
    double GetIntermediateTimeInit () ;
    void Init();

    // Results file related subs

    void WriteEndLine();
    void WriteStartLine();
    void WriteTimeToFile (double time_to_add);
    void WriteTimeToFile ();


};
#endif

#endif
